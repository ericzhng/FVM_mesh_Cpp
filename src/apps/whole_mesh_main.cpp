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
#include "polymesh/mesh_quality.hpp"
#include "polymesh/mesh_partition_manager.hpp"
#include "polymesh/local_mesh.hpp"
#include "vtkio/vtk_writer.hpp"
#include "CLI/CLI.hpp"

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

int createGeometryFromConfig(fvm::Geometry &geom, const fvm::GeometryConfig &config, double meshSize);
void assignBoundariesToEdges(const std::vector<fvm::BoundaryConfig> &boundaries,
                             const std::vector<int> &surfaceTags, bool verbose);
fvm::MeshInfo localMeshToMeshInfo(const fvm::LocalMesh &localMesh);
void writePartitionMetadata(const fvm::LocalMesh &localMesh, const std::string &fileName);
void exportMesh(const fvm::MeshInfo &meshData, const fvm::OutputConfig &output, const std::string &basePath);

// =============================================================================
// Main Function
// =============================================================================

int main(int argc, char *argv[])
{
    // 1. Setup the variables (Initialize defaults)
    std::string inputFile;
    std::string overrideOutputDir;
    bool validateOnly = false;
    bool showGui = false;
    bool verbose = false;

    // 1. HEADER: The string in the constructor becomes the description at the top.
    CLI::App app{
        "Finite Volume Mesh Generator (v1.0)"};

    // 2. FOOTER: Add examples, copyright, or contact info at the bottom.
    app.footer(
        "\nExamples:\n"
        "  ./solver channel.yaml\n"
        "  ./solver channel.yaml --output ./results --verbose\n"
        "\nCopyright (c) 2026 Eric Zhang. Distributed under MIT License.");

    // 3. Define the Flags (Booleans)
    // CLI11 automatically handles "-v" and "--verbose" mapping
    app.add_flag("--validate", validateOnly, "Validate mesh only, do not solve");
    app.add_flag("--gui", showGui, "Launch the Graphical User Interface");
    app.add_flag("-v,--verbose", verbose, "Enable verbose logging");

    // 4. Define Options (Key-Value pairs)
    app.add_option("-o,--output", overrideOutputDir, "Override the output directory");

    // 5. Define Positional Arguments (Input File)
    // By not adding a "-", this becomes a positional argument.
    // ->required() replaces your manual "if (inputFile.empty())" check.
    // ->check(CLI::ExistingFile) is a bonus: it verifies the file exists on disk!
    app.add_option("input_file", inputFile, "Path to the input mesh file")
        ->required()
        ->check(CLI::ExistingFile);

    // 6. Parse
    // This macro handles "try-catch", prints errors, prints help, and exits if needed.
    CLI11_PARSE(app, argc, argv);

    // --- Your Application Logic Starts Here ---
    std::cout << "Processing: " << inputFile << "\n";
    if (verbose)
        std::cout << "Verbose mode enabled.\n";

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
        std::cout << "   Mesh generated: " << mesher.getMeshData().elements.size() << " elements, "
                  << mesher.getMeshData().nodes.size() << " nodes\n";

        // =================================================================
        // 5. Export base mesh
        // =================================================================
        std::cout << "\n5. Exporting mesh...\n";
        const fvm::MeshInfo &meshData = mesher.getMeshData();
        std::string basePath = config.output.directory + "/" + config.output.baseName;
        exportMesh(meshData, config.output, basePath);

        // =================================================================
        // 6. Analyze mesh and write quality report
        // =================================================================
        std::string fullMshPath = config.output.directory + "/" + mshFileName;

        if (config.output.writeQualityReport || config.partition.enabled)
        {
            std::cout << "\n6. Analyzing mesh...\n";

            fvm::PolyMesh globalMesh = fvm::PolyMesh::fromGmsh(fullMshPath);
            globalMesh.analyzeMesh();

            if (verbose)
            {
                globalMesh.printSummary();
            }

            // Write markdown quality report
            if (config.output.writeQualityReport)
            {
                auto quality = fvm::MeshQuality::fromMesh(globalMesh);
                std::string qualityFile = basePath + "_quality.md";
                quality.writeMarkdownReport(qualityFile, globalMesh);
                std::cout << "   Exported: " << qualityFile << "\n";
            }

            // Partition mesh (if enabled)
            if (config.partition.enabled)
            {
                std::cout << "\n7. Partitioning mesh into " << config.partition.numParts << " parts...\n";

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
                std::cout << "\n8. Exporting partitioned meshes...\n";
                for (const auto &localMesh : localMeshes)
                {
                    std::string partitionBase = config.output.directory + "/partition_" + std::to_string(localMesh.rank);

                    // Write VTU file for visualization
                    std::string vtuFile = partitionBase + ".vtu";
                    fvm::MeshInfo partMeshInfo = localMeshToMeshInfo(localMesh);
                    fvm::VTKWriter::writeVTU(partMeshInfo, vtuFile);

                    // Write JSON metadata for parallel communication
                    if (config.output.writePartitionMetadata)
                    {
                        std::string jsonFile = partitionBase + ".json";
                        writePartitionMetadata(localMesh, jsonFile);
                    }

                    std::cout << "   Partition " << localMesh.rank << " exported\n";
                }
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

fvm::MeshInfo localMeshToMeshInfo(const fvm::LocalMesh &localMesh)
{
    fvm::MeshInfo meshInfo;

    // Nodes
    meshInfo.nodes.resize(localMesh.nNodes);
    for (std::size_t i = 0; i < localMesh.nNodes; ++i)
    {
        meshInfo.nodes[i] = localMesh.nodeCoords[i];
    }

    // Elements
    meshInfo.elements = localMesh.cellNodeConnectivity;

    // Element Types
    meshInfo.elementTypes.assign(localMesh.cellElementTypes.begin(), localMesh.cellElementTypes.end());

    return meshInfo;
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

void exportMesh(const fvm::MeshInfo &meshData, const fvm::OutputConfig &output, const std::string &basePath)
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
