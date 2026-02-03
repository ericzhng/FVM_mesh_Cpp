/**
 * @file partitioner_main.cpp
 * @brief Reads a Gmsh mesh file, partitions it, and analyzes the result.
 */

#include "polymesh/poly_mesh.hpp"
#include "polymesh/mesh_partition_manager.hpp"
#include "polymesh/local_mesh.hpp"
#include "polymesh/vtk_writer.hpp"
#include "meshgen/mesh_types.hpp"


#include <iostream>
#include <string>
#include <vector>
#include <fstream>

// Forward declaration
fvm::MeshData localMeshToMeshData(const fvm::LocalMesh& localMesh);

void printLocalMeshInfo(const fvm::LocalMesh& localMesh) {
    std::cout << "  Partition " << localMesh.rank << ":\n";
    std::cout << "    Owned Cells: " << localMesh.numOwnedCells << "\n";
    std::cout << "    Halo Cells: " << localMesh.numHaloCells << "\n";
    std::cout << "    Total Cells: " << localMesh.nCells << "\n";
    std::cout << "    Total Nodes: " << localMesh.nNodes << "\n";
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <path_to_msh_file> [num_partitions] [output_dir]\n";
        return 1;
    }

    std::string mshFile = argv[1];
    int nParts = 4;
    if (argc > 2) {
        nParts = std::stoi(argv[2]);
    }
    std::string outputDir = "data";
    if (argc > 3) {
        outputDir = argv[3];
    }


    std::cout << "FVM Mesh Partitioner\n";
    std::cout << "====================\n";
    std::cout << "Parameters:\n";
    std::cout << "  Mesh file: " << mshFile << "\n";
    std::cout << "  Partitions: " << nParts << "\n";
    std::cout << "  Output dir: " << outputDir << "\n\n";

    // Check if the mesh file exists
    std::ifstream f(mshFile.c_str());
    if (!f.good()) {
        std::cerr << "Error: Mesh file not found or not readable: " << mshFile << std::endl;
        return 1;
    }
    f.close();

    try {
        // 1. Create a PolyMesh from the Gmsh file
        std::cout << "1. Reading and analyzing global mesh...\n";
        fvm::PolyMesh globalMesh = fvm::PolyMesh::fromGmsh(mshFile);
        globalMesh.analyzeMesh();
        globalMesh.printSummary();

        // 2. Partition the mesh
        std::cout << "\n2. Partitioning mesh into " << nParts << " domains...\n";
        std::vector<fvm::LocalMesh> localMeshes = fvm::MeshPartitionManager::createLocalMeshes(globalMesh, nParts, "metis");
        std::cout << "Partitioning complete.\n\n";

        // 3. Analyze the local meshes
        std::cout << "3. Analyzing local meshes:\n";
        for (const auto& localMesh : localMeshes) {
            printLocalMeshInfo(localMesh);
        }

        // 4. Write local meshes to VTU files
        std::cout << "\n4. Writing local meshes to VTU files...\n";
        for (const auto& localMesh : localMeshes) {
            std::string fileName = outputDir + "/partition_" + std::to_string(localMesh.rank) + ".vtu";
            fvm::MeshData meshData = localMeshToMeshData(localMesh);
            fvm::VTKWriter::writeVTU(meshData, fileName);
            std::cout << "  - Wrote " << fileName << "\n";
        }

        std::cout << "\nDone!\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
}

fvm::MeshData localMeshToMeshData(const fvm::LocalMesh& localMesh) {
    fvm::MeshData meshData;

    // Nodes
    meshData.nodes.resize(localMesh.nNodes);
    meshData.nodeIds.resize(localMesh.nNodes);
    for (size_t i = 0; i < localMesh.nNodes; ++i) {
        meshData.nodes[i] = localMesh.nodeCoords[i];
        meshData.nodeIds[i] = localMesh.l2gNodes[i];
    }

    // Cells
    meshData.cells = localMesh.cellNodeConnectivity;

    // Cell Types
    for (int elemType : localMesh.cellElementTypes) {
        auto it = localMesh.elementTypeProperties.find(elemType);
        if (it != localMesh.elementTypeProperties.end()) {
            meshData.cellTypes.push_back(fvm::getVTKCellType(it->second.numNodes));
        } else {
            // Fallback for unknown types
            meshData.cellTypes.push_back(fvm::VTKCellType::POLYGON);
        }
    }
    
    // Boundary information can be tricky to reconstruct perfectly without more
    // context from the global mesh, but we can approximate it for visualization.
    // This implementation will be simplified.

    return meshData;
}
