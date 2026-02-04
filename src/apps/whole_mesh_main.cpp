/**
 * @file whole_mesh_main.cpp
 * @brief Unified mesh generation executable with YAML input.
 *
 * This executable reads a YAML configuration file and performs:
 * - Geometry creation (multiple shapes with union)
 * - Mesh generation
 * - Boundary identification via expressions
 * - Optional mesh partitioning with reordering
 * - Export to multiple formats
 */

#include "common/fvm_types.hpp"
#include "input/input_config.hpp"
#include "input/input_parser.hpp"
#include "input/expr_evaluator.hpp"
#include "meshgen/geometry.hpp"
#include "meshgen/mesh_generator.hpp"
#include "polymesh/poly_mesh.hpp"
#include "polymesh/mesh_partition_manager.hpp"
#include "polymesh/local_mesh.hpp"
#include "vtkio/vtk_writer.hpp"

#include <gmsh.h>

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>

namespace fs = std::filesystem;

// =============================================================================
// Forward Declarations
// =============================================================================

void printUsage(const char *progName);
int createGeometryFromConfig(fvm::Geometry &geom, const fvm::GeometryConfig &config, double meshSize);
void assignBoundariesToEdges(const std::vector<fvm::BoundaryConfig> &boundaries,
                             const std::vector<int> &surfaceTags, bool verbose);
fvm::MeshData localMeshToMeshData(const fvm::LocalMesh &localMesh);
void writePartitionMetadata(const fvm::LocalMesh &localMesh, const std::string &fileName);
void exportMesh(const fvm::MeshData &meshData, const fvm::OutputConfig &output, const std::string &basePath);

// =============================================================================
// Main Function
// =============================================================================

int main(int argc, char *argv[])
{
    // Parse command-line arguments
    std::string inputFile;
    std::string overrideOutputDir;
    bool validateOnly = false;
    bool showGui = false;
    bool verbose = false;

    for (int i = 1; i < argc; ++i)
    {
        std::string arg = argv[i];

        if (arg == "--help" || arg == "-h")
        {
            printUsage(argv[0]);
            return 0;
        }
        else if (arg == "--validate")
        {
            validateOnly = true;
        }
        else if (arg == "--gui")
        {
            showGui = true;
        }
        else if (arg == "--verbose" || arg == "-v")
        {
            verbose = true;
        }
        else if ((arg == "--output" || arg == "-o") && i + 1 < argc)
        {
            overrideOutputDir = argv[++i];
        }
        else if (arg[0] != '-')
        {
            inputFile = arg;
        }
        else
        {
            std::cerr << "Unknown option: " << arg << "\n";
            printUsage(argv[0]);
            return 1;
        }
    }

    if (inputFile.empty())
    {
        std::cerr << "Error: No input file specified.\n";
        printUsage(argv[0]);
        return 1;
    }

    // Check if input file exists
    if (!fs::exists(inputFile))
    {
        std::cerr << "Error: Input file not found: " << inputFile << "\n";
        return 1;
    }

    std::cout << "FVM Unified Mesh Generator\n";
    std::cout << "==========================\n";
    std::cout << "Input file: " << inputFile << "\n\n";

    try
    {
        // =================================================================
        // 1. Parse and validate configuration
        // =================================================================
        std::cout << "1. Parsing configuration...\n";
        fvm::InputConfig config = fvm::parseYamlFile(inputFile);

        // Override output directory if specified
        if (!overrideOutputDir.empty())
        {
            config.output.directory = overrideOutputDir;
        }

        // Validate configuration
        std::string validationError;
        if (!config.validate(validationError))
        {
            std::cerr << "Configuration validation failed:\n"
                      << validationError;
            return 1;
        }
        std::cout << "   Configuration valid.\n";

        if (verbose)
        {
            std::cout << "   Project: " << config.projectName << "\n";
            std::cout << "   Geometries: " << config.geometries.size() << "\n";
            std::cout << "   Boundaries: " << config.boundaries.size() << "\n";
            std::cout << "   Partition: " << (config.partition.enabled ? "enabled" : "disabled") << "\n";
        }

        if (validateOnly)
        {
            std::cout << "\nValidation complete. Exiting.\n";
            return 0;
        }

        // Create output directory
        fs::path outputPath(config.output.directory);
        if (!fs::exists(outputPath))
        {
            fs::create_directories(outputPath);
        }

        // =================================================================
        // 2. Initialize Gmsh and create geometries
        // =================================================================
        std::cout << "\n2. Creating geometry...\n";
        gmsh::initialize();
        gmsh::model::add(config.projectName);

        fvm::Geometry geom(config.projectName);
        std::vector<int> surfaceTags;

        for (std::size_t i = 0; i < config.geometries.size(); ++i)
        {
            const auto &geoConfig = config.geometries[i];
            int tag = createGeometryFromConfig(geom, geoConfig, config.mesh.meshSize);

            if (tag <= 0)
            {
                std::cerr << "Error: Failed to create geometry " << i << " (" << geoConfig.type << ")\n";
                gmsh::finalize();
                return 1;
            }

            surfaceTags.push_back(tag);
            std::cout << "   Created " << geoConfig.type;
            if (!geoConfig.id.empty())
            {
                std::cout << " '" << geoConfig.id << "'";
            }
            std::cout << " (surface tag " << tag << ")\n";
        }

        // Synchronize the geometry kernel
        gmsh::model::geo::synchronize();

        // =================================================================
        // 3. Assign boundaries using expressions
        // =================================================================
        std::cout << "\n3. Assigning boundaries...\n";
        assignBoundariesToEdges(config.boundaries, surfaceTags, verbose);

        // =================================================================
        // 4. Generate mesh
        // =================================================================
        std::cout << "\n4. Generating mesh...\n";
        fvm::MeshGenerator mesher(surfaceTags, config.output.directory);

        std::map<int, fvm::MeshParams> meshParams;
        for (int tag : surfaceTags)
        {
            meshParams[tag] = {config.mesh.meshType, config.mesh.charLength};
        }

        std::string mshFileName = config.output.baseName + ".msh";
        mesher.generate(meshParams, mshFileName);
        std::cout << "   Mesh generated: " << mesher.getMeshData().cells.size() << " cells, "
                  << mesher.getMeshData().nodes.size() << " nodes\n";

        // =================================================================
        // 5. Export base mesh
        // =================================================================
        std::cout << "\n5. Exporting mesh...\n";
        const fvm::MeshData &meshData = mesher.getMeshData();
        std::string basePath = config.output.directory + "/" + config.output.baseName;
        exportMesh(meshData, config.output, basePath);

        // =================================================================
        // 6. Partition mesh (if enabled)
        // =================================================================
        if (config.partition.enabled)
        {
            std::cout << "\n6. Partitioning mesh into " << config.partition.numParts << " parts...\n";

            // Read mesh into PolyMesh
            std::string fullMshPath = config.output.directory + "/" + mshFileName;
            fvm::PolyMesh globalMesh = fvm::PolyMesh::fromGmsh(fullMshPath);
            globalMesh.analyzeMesh();

            if (verbose)
            {
                globalMesh.printSummary();
            }

            // Partition using MeshPartitionManager
            std::vector<fvm::LocalMesh> localMeshes = fvm::MeshPartitionManager::createLocalMeshes(
                globalMesh,
                config.partition.numParts,
                config.partition.method,
                config.reorder.cellStrategy,
                config.reorder.nodeStrategy);

            std::cout << "   Partitioning complete.\n";

            // Analyze partitions
            for (const auto &localMesh : localMeshes)
            {
                std::cout << "   Partition " << localMesh.rank << ": "
                          << localMesh.numOwnedCells << " owned, "
                          << localMesh.numHaloCells << " halo cells\n";
            }

            // Export partitioned meshes
            std::cout << "\n7. Exporting partitioned meshes...\n";
            for (const auto &localMesh : localMeshes)
            {
                std::string partitionBase = config.output.directory + "/partition_" + std::to_string(localMesh.rank);

                // Write VTU file for visualization
                std::string vtuFile = partitionBase + ".vtu";
                fvm::MeshData partMeshData = localMeshToMeshData(localMesh);
                fvm::VTKWriter::writeVTU(partMeshData, vtuFile);

                // Write JSON metadata for parallel communication
                if (config.output.writePartitionMetadata)
                {
                    std::string jsonFile = partitionBase + ".json";
                    writePartitionMetadata(localMesh, jsonFile);
                }

                std::cout << "   Partition " << localMesh.rank << " exported\n";
            }
        }

        // Show GUI if requested
        if (showGui)
        {
            std::cout << "\nOpening Gmsh GUI...\n";
            gmsh::fltk::run();
        }

        // Cleanup
        gmsh::finalize();

        std::cout << "\nDone!\n";
        return 0;
    }
    catch (const std::exception &e)
    {
        std::cerr << "Error: " << e.what() << std::endl;
        try
        {
            gmsh::finalize();
        }
        catch (...)
        {
        }
        return 1;
    }
}

// =============================================================================
// Helper Functions
// =============================================================================

void printUsage(const char *progName)
{
    std::cout << "Usage: " << progName << " <input.yaml> [options]\n\n"
              << "Options:\n"
              << "  --output, -o <dir>    Override output directory\n"
              << "  --validate            Validate input file only\n"
              << "  --gui                 Show Gmsh GUI after generation\n"
              << "  --verbose, -v         Enable verbose output\n"
              << "  --help, -h            Show this help message\n\n"
              << "Example:\n"
              << "  " << progName << " channel.yaml\n"
              << "  " << progName << " channel.yaml --output ./results --verbose\n";
}

int createGeometryFromConfig(fvm::Geometry &geom, const fvm::GeometryConfig &config, double meshSize)
{
    if (config.type == "rectangle")
    {
        return geom.rectangle(config.length, config.width, config.x, config.y, meshSize);
    }
    else if (config.type == "circle")
    {
        return geom.circle(config.radius, config.x, config.y, meshSize);
    }
    else if (config.type == "triangle")
    {
        return geom.triangle(config.p1, config.p2, config.p3, meshSize);
    }
    else if (config.type == "ellipse")
    {
        return geom.ellipse(config.r1, config.r2, config.x, config.y, meshSize);
    }
    else if (config.type == "polygon")
    {
        return geom.polygon(config.points, config.convexHull, meshSize);
    }

    return -1; // Unknown type
}

void assignBoundariesToEdges(const std::vector<fvm::BoundaryConfig> &boundaries,
                             const std::vector<int> &surfaceTags, bool verbose)
{
    // Get all boundary curves from all surfaces
    std::vector<std::pair<int, int>> allEdges;
    for (int surfTag : surfaceTags)
    {
        std::vector<std::pair<int, int>> edges;
        gmsh::model::getBoundary({{2, surfTag}}, edges);
        allEdges.insert(allEdges.end(), edges.begin(), edges.end());
    }

    if (verbose)
    {
        std::cout << "   Found " << allEdges.size() << " boundary edges\n";
    }

    // Process each boundary definition
    for (const auto &bc : boundaries)
    {
        fvm::ExprEvaluator eval(bc.expr);

        if (!eval.isValid())
        {
            std::cerr << "   Warning: Invalid expression for boundary '" << bc.name
                      << "': " << eval.getError() << "\n";
            continue;
        }

        std::vector<int> matchingEdges;

        for (const auto &[dim, tag] : allEdges)
        {
            // Get edge bounding box to compute midpoint
            double minX, minY, minZ, maxX, maxY, maxZ;
            gmsh::model::getBoundingBox(1, std::abs(tag), minX, minY, minZ, maxX, maxY, maxZ);

            double midX = (minX + maxX) / 2.0;
            double midY = (minY + maxY) / 2.0;

            // Evaluate expression at midpoint
            if (eval.matches(midX, midY))
            {
                matchingEdges.push_back(std::abs(tag));
            }
        }

        if (!matchingEdges.empty())
        {
            gmsh::model::addPhysicalGroup(1, matchingEdges, -1, bc.name);
            std::cout << "   Boundary '" << bc.name << "': " << matchingEdges.size() << " edges\n";
        }
        else
        {
            std::cout << "   Warning: No edges matched for boundary '" << bc.name << "'\n";
        }
    }

    // Create a physical group for the surfaces
    gmsh::model::addPhysicalGroup(2, surfaceTags, -1, "domain");
}

fvm::MeshData localMeshToMeshData(const fvm::LocalMesh &localMesh)
{
    fvm::MeshData meshData;

    // Nodes
    meshData.nodes.resize(localMesh.nNodes);
    meshData.nodeIds.resize(localMesh.nNodes);
    for (std::size_t i = 0; i < localMesh.nNodes; ++i)
    {
        meshData.nodes[i] = localMesh.nodeCoords[i];
        meshData.nodeIds[i] = localMesh.l2gNodes[i];
    }

    // Cells
    meshData.cells = localMesh.cellNodeConnectivity;

    // Cell Types
    for (int elemType : localMesh.cellElementTypes)
    {
        auto it = localMesh.elementTypeProperties.find(elemType);
        if (it != localMesh.elementTypeProperties.end())
        {
            meshData.cellTypes.push_back(fvm::getVTKCellType(it->second.numNodes));
        }
        else
        {
            meshData.cellTypes.push_back(fvm::VTKCellType::POLYGON);
        }
    }

    return meshData;
}

void writePartitionMetadata(const fvm::LocalMesh &localMesh, const std::string &fileName)
{
    std::ofstream ofs(fileName);
    if (!ofs.is_open())
    {
        throw std::runtime_error("Failed to open file for writing: " + fileName);
    }

    ofs << "{\n";

    // Basic partition info
    ofs << "  \"rank\": " << localMesh.rank << ",\n";
    ofs << "  \"numOwnedCells\": " << localMesh.numOwnedCells << ",\n";
    ofs << "  \"numHaloCells\": " << localMesh.numHaloCells << ",\n";
    ofs << "  \"totalCells\": " << localMesh.nCells << ",\n";
    ofs << "  \"totalNodes\": " << localMesh.nNodes << ",\n";

    // Local-to-global cell mapping
    ofs << "  \"l2gCells\": [";
    for (std::size_t i = 0; i < localMesh.l2gCells.size(); ++i)
    {
        if (i > 0)
            ofs << ", ";
        ofs << localMesh.l2gCells[i];
    }
    ofs << "],\n";

    // Local-to-global node mapping
    ofs << "  \"l2gNodes\": [";
    for (std::size_t i = 0; i < localMesh.l2gNodes.size(); ++i)
    {
        if (i > 0)
            ofs << ", ";
        ofs << localMesh.l2gNodes[i];
    }
    ofs << "],\n";

    // Send map
    ofs << "  \"sendMap\": {\n";
    {
        bool firstRank = true;
        for (const auto &[neighborRank, cellIndices] : localMesh.sendMap)
        {
            if (!firstRank)
                ofs << ",\n";
            firstRank = false;
            ofs << "    \"" << neighborRank << "\": [";
            for (std::size_t i = 0; i < cellIndices.size(); ++i)
            {
                if (i > 0)
                    ofs << ", ";
                ofs << cellIndices[i];
            }
            ofs << "]";
        }
        ofs << "\n  },\n";
    }

    // Receive map
    ofs << "  \"recvMap\": {\n";
    {
        bool firstRank = true;
        for (const auto &[neighborRank, cellIndices] : localMesh.recvMap)
        {
            if (!firstRank)
                ofs << ",\n";
            firstRank = false;
            ofs << "    \"" << neighborRank << "\": [";
            for (std::size_t i = 0; i < cellIndices.size(); ++i)
            {
                if (i > 0)
                    ofs << ", ";
                ofs << cellIndices[i];
            }
            ofs << "]";
        }
        ofs << "\n  }\n";
    }

    ofs << "}\n";
    ofs.close();
}

void exportMesh(const fvm::MeshData &meshData, const fvm::OutputConfig &output, const std::string &basePath)
{
    for (const auto &format : output.formats)
    {
        if (format == "vtk")
        {
            std::string fileName = basePath + ".vtk";
            fvm::VTKWriter::writeVTK(meshData, fileName);
            std::cout << "   Exported: " << fileName << "\n";
        }
        else if (format == "vtu")
        {
            std::string fileName = basePath + ".vtu";
            fvm::VTKWriter::writeVTU(meshData, fileName);
            std::cout << "   Exported: " << fileName << "\n";
        }
        else if (format == "openfoam")
        {
            std::string dirName = basePath + "_openfoam";
            fvm::VTKWriter::writeOpenFOAM(meshData, dirName);
            std::cout << "   Exported: " << dirName << "/\n";
        }
        // msh is already written by MeshGenerator
    }

    if (output.writeBoundaryInfo)
    {
        std::string fileName = basePath + "_boundaries.txt";
        fvm::VTKWriter::writeBoundaryInfo(meshData, fileName);
        std::cout << "   Exported: " << fileName << "\n";
    }
}
