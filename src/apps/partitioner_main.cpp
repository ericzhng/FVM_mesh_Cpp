/**
 * @file partitioner_main.cpp
 * @brief Reads a Gmsh mesh file, partitions it, and analyzes the result.
 */

#include "common/fvm_types.hpp"
#include "vtkio/vtk_writer.hpp"
#include "polymesh/poly_mesh.hpp"
#include "polymesh/mesh_partition_manager.hpp"
#include "polymesh/local_mesh.hpp"
#include <gmsh.h>

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <filesystem>

// Forward declarations
fvm::MeshInfo localMeshToMeshInfo(const fvm::LocalMesh &localMesh);
void writePartitionMetadata(const fvm::LocalMesh &localMesh, const std::string &fileName);

void printLocalMeshInfo(const fvm::LocalMesh &localMesh)
{
    std::cout << "  Partition " << localMesh.rank << ":\n";
    std::cout << "    Owned Cells: " << localMesh.numOwnedCells << "\n";
    std::cout << "    Halo Cells: " << localMesh.numHaloCells << "\n";
    std::cout << "    Total Cells: " << localMesh.nCells << "\n";
    std::cout << "    Total Nodes: " << localMesh.nNodes << "\n";
}

int main(int argc, char *argv[])
{
    if (argc < 2)
    {
        std::cerr << "Usage: " << argv[0] << " <path_to_msh_file> [num_partitions] [output_dir] [reorder_strategy]\n";
        std::cerr << "  reorder_strategy: rcm (default), gps, sloan, spectral, spatial_x, spatial_y, random, none\n";
        return 1;
    }

    std::string mshFile = argv[1];
    int nParts = 4;
    if (argc > 2)
    {
        nParts = std::stoi(argv[2]);
    }
    std::string outputDir = "data";
    if (argc > 3)
    {
        outputDir = argv[3];
    }
    std::string reorderStrategy = "rcm";
    if (argc > 4)
    {
        reorderStrategy = argv[4];
    }

    std::filesystem::path outputPath(outputDir);
    if (!std::filesystem::exists(outputPath))
    {
        std::filesystem::create_directories(outputPath);
    }

    std::cout << "FVM Mesh Partitioner\n";
    std::cout << "====================\n";
    std::cout << "Parameters:\n";
    std::cout << "  Mesh file: " << mshFile << "\n";
    std::cout << "  Partitions: " << nParts << "\n";
    std::cout << "  Output dir: " << outputDir << "\n";
    std::cout << "  Reorder strategy: " << reorderStrategy << "\n\n";

    // Check if the mesh file exists
    std::ifstream f(mshFile.c_str());
    if (!f.good())
    {
        std::cerr << "Error: Mesh file not found or not readable: " << mshFile << std::endl;
        return 1;
    }
    f.close();

    gmsh::initialize();

    try
    {
        // 1. Create a PolyMesh from the Gmsh file
        std::cout << "1. Reading and analyzing global mesh...\n";
        fvm::PolyMesh globalMesh = fvm::PolyMesh::fromGmsh(mshFile);
        globalMesh.analyzeMesh();
        globalMesh.printSummary();

        // 2. Partition the mesh
        std::cout << "\n2. Partitioning mesh into " << nParts << " domains...\n";
        std::vector<fvm::LocalMesh> localMeshes = fvm::MeshPartitionManager::createLocalMeshes(globalMesh, nParts, "metis");
        std::cout << "Partitioning complete.\n\n";

        // 3. Reorder cells and nodes for each partition
        if (reorderStrategy != "none")
        {
            std::cout << "3. Reordering cells and nodes (strategy: " << reorderStrategy << ")...\n";
            for (auto &localMesh : localMeshes)
            {
                localMesh.reorderCells(reorderStrategy);
                localMesh.reorderNodes(reorderStrategy);
                std::cout << "  - Partition " << localMesh.rank << " reordered\n";
            }
            std::cout << "Reordering complete.\n\n";
        }
        else
        {
            std::cout << "3. Skipping reordering (strategy: none)\n\n";
        }

        // 4. Analyze the local meshes
        std::cout << "4. Analyzing local meshes:\n";
        for (const auto &localMesh : localMeshes)
        {
            printLocalMeshInfo(localMesh);
        }

        // 5. Write local meshes to VTU and metadata to JSON files
        std::cout << "\n5. Writing local meshes to VTU and metadata to JSON files...\n";
        for (const auto &localMesh : localMeshes)
        {
            std::string baseName = outputDir + "/partition_" + std::to_string(localMesh.rank);

            // Write VTU file for visualization
            std::string vtuFileName = baseName + ".vtu";
            fvm::MeshInfo meshData = localMeshToMeshInfo(localMesh);
            fvm::VTKWriter::writeVTU(meshData, vtuFileName);

            // Write JSON metadata for parallel communication
            std::string jsonFileName = baseName + ".json";
            writePartitionMetadata(localMesh, jsonFileName);

            std::cout << "  - Partition " << localMesh.rank << ": " << vtuFileName << ", " << jsonFileName << "\n";
        }

        std::cout << "\nDone!\n";
        gmsh::finalize();
        return 0;
    }
    catch (const std::exception &e)
    {
        std::cerr << "Error: " << e.what() << std::endl;
        gmsh::finalize();
        return 1;
    }
}

fvm::MeshInfo localMeshToMeshInfo(const fvm::LocalMesh &localMesh)
{
    fvm::MeshInfo meshInfo;

    // Nodes
    meshInfo.nodes.resize(localMesh.nNodes);
    for (size_t i = 0; i < localMesh.nNodes; ++i)
    {
        meshInfo.nodes[i] = localMesh.nodeCoords[i];
    }

    // Elements
    meshInfo.elements = localMesh.cellNodeConnectivity;

    // Element Types
    for (int elemType : localMesh.cellElementTypes)
    {
        auto it = localMesh.elementTypeProperties.find(elemType);
        if (it != localMesh.elementTypeProperties.end())
        {
            meshInfo.elementTypes.push_back(fvm::getVTKCellType(it->second.numNodes));
        }
        else
        {
            // Fallback for unknown types
            meshInfo.elementTypes.push_back(fvm::VTKCellType::POLYGON);
        }
    }

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
    for (size_t i = 0; i < localMesh.l2gCells.size(); ++i)
    {
        if (i > 0) ofs << ", ";
        ofs << localMesh.l2gCells[i];
    }
    ofs << "],\n";

    // Local-to-global node mapping
    ofs << "  \"l2gNodes\": [";
    for (size_t i = 0; i < localMesh.l2gNodes.size(); ++i)
    {
        if (i > 0) ofs << ", ";
        ofs << localMesh.l2gNodes[i];
    }
    ofs << "],\n";

    // Send map: cells to send to each neighbor rank
    ofs << "  \"sendMap\": {\n";
    {
        bool firstRank = true;
        for (const auto &[neighborRank, cellIndices] : localMesh.sendMap)
        {
            if (!firstRank) ofs << ",\n";
            firstRank = false;
            ofs << "    \"" << neighborRank << "\": [";
            for (size_t i = 0; i < cellIndices.size(); ++i)
            {
                if (i > 0) ofs << ", ";
                ofs << cellIndices[i];
            }
            ofs << "]";
        }
        ofs << "\n  },\n";
    }

    // Receive map: halo cells to receive from each neighbor rank
    ofs << "  \"recvMap\": {\n";
    {
        bool firstRank = true;
        for (const auto &[neighborRank, cellIndices] : localMesh.recvMap)
        {
            if (!firstRank) ofs << ",\n";
            firstRank = false;
            ofs << "    \"" << neighborRank << "\": [";
            for (size_t i = 0; i < cellIndices.size(); ++i)
            {
                if (i > 0) ofs << ", ";
                ofs << cellIndices[i];
            }
            ofs << "]";
        }
        ofs << "\n  }\n";
    }

    ofs << "}\n";
    ofs.close();
}
