#include <gtest/gtest.h>
#include "polymesh/poly_mesh.hpp"
#include "polymesh/mesh_quality.hpp"
#include "polymesh/partition.hpp"
#include <cmath>

namespace fvm
{
    namespace testing
    {

        // =============================================================================
        // PolyMesh Factory Method Tests
        // =============================================================================

        TEST(PolyMeshFactory, CreateStructuredQuadMesh2x2)
        {
            auto mesh = PolyMesh::createStructuredQuadMesh(2, 2);

            EXPECT_EQ(mesh.dimension, 2);
            EXPECT_EQ(mesh.nCells, 4u);
            EXPECT_EQ(mesh.nNodes, 9u); // (2+1) * (2+1)
            EXPECT_TRUE(mesh.isAnalyzed());
        }

        TEST(PolyMeshFactory, CreateStructuredQuadMesh3x4)
        {
            auto mesh = PolyMesh::createStructuredQuadMesh(3, 4);

            EXPECT_EQ(mesh.dimension, 2);
            EXPECT_EQ(mesh.nCells, 12u); // 3 * 4
            EXPECT_EQ(mesh.nNodes, 20u); // (3+1) * (4+1)
        }

        TEST(PolyMeshFactory, CreateStructuredQuadMesh1x1)
        {
            auto mesh = PolyMesh::createStructuredQuadMesh(1, 1);

            EXPECT_EQ(mesh.nCells, 1u);
            EXPECT_EQ(mesh.nNodes, 4u);
        }

        // =============================================================================
        // Node Coordinate Tests
        // =============================================================================

        TEST(PolyMeshNodes, NodeCoordinatesSize)
        {
            auto mesh = PolyMesh::createStructuredQuadMesh(3, 3);

            EXPECT_EQ(mesh.nodeCoords.size(), mesh.nNodes);
        }

        TEST(PolyMeshNodes, NodeCoordinatesBoundingBox)
        {
            auto mesh = PolyMesh::createStructuredQuadMesh(3, 3);

            // For a unit square mesh
            double minX = std::numeric_limits<double>::max();
            double maxX = std::numeric_limits<double>::lowest();
            double minY = std::numeric_limits<double>::max();
            double maxY = std::numeric_limits<double>::lowest();

            for (const auto &coord : mesh.nodeCoords)
            {
                minX = std::min(minX, coord[0]);
                maxX = std::max(maxX, coord[0]);
                minY = std::min(minY, coord[1]);
                maxY = std::max(maxY, coord[1]);
            }

            EXPECT_NEAR(minX, 0.0, 1e-10);
            EXPECT_NEAR(maxX, 1.0, 1e-10);
            EXPECT_NEAR(minY, 0.0, 1e-10);
            EXPECT_NEAR(maxY, 1.0, 1e-10);
        }

        // =============================================================================
        // Cell Connectivity Tests
        // =============================================================================

        TEST(PolyMeshCells, CellConnectivitySize)
        {
            auto mesh = PolyMesh::createStructuredQuadMesh(3, 3);

            EXPECT_EQ(mesh.cellNodeConnectivity.size(), mesh.nCells);
        }

        TEST(PolyMeshCells, QuadCellsHaveFourNodes)
        {
            auto mesh = PolyMesh::createStructuredQuadMesh(3, 3);

            for (const auto &cell : mesh.cellNodeConnectivity)
            {
                EXPECT_EQ(cell.size(), 4u);
            }
        }

        TEST(PolyMeshCells, CellNodeIndicesValid)
        {
            auto mesh = PolyMesh::createStructuredQuadMesh(3, 3);

            for (const auto &cell : mesh.cellNodeConnectivity)
            {
                for (auto nodeIdx : cell)
                {
                    EXPECT_LT(nodeIdx, mesh.nNodes);
                }
            }
        }

        // =============================================================================
        // Cell Centroid Tests
        // =============================================================================

        TEST(PolyMeshCentroids, CentroidsSize)
        {
            auto mesh = PolyMesh::createStructuredQuadMesh(3, 3);

            EXPECT_EQ(mesh.cellCentroids.size(), mesh.nCells);
        }

        TEST(PolyMeshCentroids, CentroidsInsideBounds)
        {
            auto mesh = PolyMesh::createStructuredQuadMesh(3, 3);

            for (const auto &centroid : mesh.cellCentroids)
            {
                EXPECT_GE(centroid[0], 0.0);
                EXPECT_LE(centroid[0], 1.0);
                EXPECT_GE(centroid[1], 0.0);
                EXPECT_LE(centroid[1], 1.0);
            }
        }

        // =============================================================================
        // Cell Volume Tests
        // =============================================================================

        TEST(PolyMeshVolumes, VolumesSize)
        {
            auto mesh = PolyMesh::createStructuredQuadMesh(3, 3);

            EXPECT_EQ(mesh.cellVolumes.size(), mesh.nCells);
        }

        TEST(PolyMeshVolumes, VolumesPositive)
        {
            auto mesh = PolyMesh::createStructuredQuadMesh(3, 3);

            for (double vol : mesh.cellVolumes)
            {
                EXPECT_GT(vol, 0.0);
            }
        }

        TEST(PolyMeshVolumes, TotalVolumeEqualsOne)
        {
            auto mesh = PolyMesh::createStructuredQuadMesh(5, 5);

            double totalVolume = 0.0;
            for (double vol : mesh.cellVolumes)
            {
                totalVolume += vol;
            }

            // Unit square mesh should have total area of 1.0
            EXPECT_NEAR(totalVolume, 1.0, 1e-10);
        }

        TEST(PolyMeshVolumes, UniformCellSizes)
        {
            auto mesh = PolyMesh::createStructuredQuadMesh(4, 4);

            double expectedVolume = 1.0 / 16.0; // 1/(4*4)
            for (double vol : mesh.cellVolumes)
            {
                EXPECT_NEAR(vol, expectedVolume, 1e-10);
            }
        }

        // =============================================================================
        // Cell Neighbor Tests
        // =============================================================================

        TEST(PolyMeshNeighbors, NeighborsSize)
        {
            auto mesh = PolyMesh::createStructuredQuadMesh(3, 3);

            EXPECT_EQ(mesh.cellNeighbors.size(), mesh.nCells);
        }

        TEST(PolyMeshNeighbors, QuadCellsHaveFourFaces)
        {
            auto mesh = PolyMesh::createStructuredQuadMesh(3, 3);

            for (const auto &neighbors : mesh.cellNeighbors)
            {
                EXPECT_EQ(neighbors.size(), 4u);
            }
        }

        TEST(PolyMeshNeighbors, InteriorCellHasFourNeighbors)
        {
            auto mesh = PolyMesh::createStructuredQuadMesh(3, 3);

            // Cell in the center (1,1) in a 3x3 grid should have 4 neighbors
            // Cell index 4 is the center cell (row 1, col 1) in row-major order
            int centerCellIdx = 4;
            int neighborCount = 0;
            for (int neighbor : mesh.cellNeighbors[centerCellIdx])
            {
                if (neighbor >= 0)
                {
                    neighborCount++;
                }
            }
            EXPECT_EQ(neighborCount, 4);
        }

        TEST(PolyMeshNeighbors, CornerCellHasTwoNeighbors)
        {
            auto mesh = PolyMesh::createStructuredQuadMesh(3, 3);

            // Bottom-left corner cell (index 0) should have only 2 neighbors
            int neighborCount = 0;
            for (int neighbor : mesh.cellNeighbors[0])
            {
                if (neighbor >= 0)
                {
                    neighborCount++;
                }
            }
            EXPECT_EQ(neighborCount, 2);
        }

        TEST(PolyMeshNeighbors, NeighborSymmetry)
        {
            auto mesh = PolyMesh::createStructuredQuadMesh(3, 3);

            // If cell A is a neighbor of cell B, then B should be a neighbor of A
            for (std::size_t cellA = 0; cellA < mesh.nCells; ++cellA)
            {
                for (int neighborB : mesh.cellNeighbors[cellA])
                {
                    if (neighborB >= 0)
                    {
                        bool found = false;
                        for (int neighborOfB : mesh.cellNeighbors[neighborB])
                        {
                            if (neighborOfB == static_cast<int>(cellA))
                            {
                                found = true;
                                break;
                            }
                        }
                        EXPECT_TRUE(found) << "Cell " << cellA << " has neighbor " << neighborB
                                           << " but " << neighborB << " doesn't have " << cellA << " as neighbor";
                    }
                }
            }
        }

        // =============================================================================
        // Face Properties Tests
        // =============================================================================

        TEST(PolyMeshFaces, FaceAreasSize)
        {
            auto mesh = PolyMesh::createStructuredQuadMesh(3, 3);

            EXPECT_EQ(mesh.cellFaceAreas.size(), mesh.nCells);
        }

        TEST(PolyMeshFaces, FaceAreasPositive)
        {
            auto mesh = PolyMesh::createStructuredQuadMesh(3, 3);

            for (const auto &faceAreas : mesh.cellFaceAreas)
            {
                for (double area : faceAreas)
                {
                    EXPECT_GT(area, 0.0);
                }
            }
        }

        TEST(PolyMeshFaces, FaceNormalsSize)
        {
            auto mesh = PolyMesh::createStructuredQuadMesh(3, 3);

            EXPECT_EQ(mesh.cellFaceNormals.size(), mesh.nCells);
        }

        TEST(PolyMeshFaces, FaceNormalsUnitLength)
        {
            auto mesh = PolyMesh::createStructuredQuadMesh(3, 3);

            for (const auto &faceNormals : mesh.cellFaceNormals)
            {
                for (const auto &normal : faceNormals)
                {
                    double length = std::sqrt(normal[0] * normal[0] +
                                              normal[1] * normal[1] +
                                              normal[2] * normal[2]);
                    EXPECT_NEAR(length, 1.0, 1e-10);
                }
            }
        }

        // =============================================================================
        // Mesh Quality Tests
        // =============================================================================

        TEST(MeshQualityTest, ComputeFromMesh)
        {
            auto mesh = PolyMesh::createStructuredQuadMesh(3, 3);
            auto quality = MeshQuality::fromMesh(mesh);

            EXPECT_GT(quality.minMaxVolumeRatio, 0.0);
            EXPECT_LE(quality.minMaxVolumeRatio, 1.0);
        }

        TEST(MeshQualityTest, SkewnessSize)
        {
            auto mesh = PolyMesh::createStructuredQuadMesh(3, 3);
            auto quality = MeshQuality::fromMesh(mesh);

            EXPECT_EQ(quality.cellSkewness.size(), mesh.nCells);
        }

        TEST(MeshQualityTest, SkewnessRange)
        {
            auto mesh = PolyMesh::createStructuredQuadMesh(3, 3);
            auto quality = MeshQuality::fromMesh(mesh);

            for (double skewness : quality.cellSkewness)
            {
                EXPECT_GE(skewness, 0.0);
                EXPECT_LE(skewness, 1.0);
            }
        }

        TEST(MeshQualityTest, UniformMeshLowSkewness)
        {
            auto mesh = PolyMesh::createStructuredQuadMesh(5, 5);
            auto quality = MeshQuality::fromMesh(mesh);

            // A uniform structured mesh should have very low skewness
            for (double skewness : quality.cellSkewness)
            {
                EXPECT_LT(skewness, 0.1);
            }
        }

        TEST(MeshQualityTest, AspectRatioSize)
        {
            auto mesh = PolyMesh::createStructuredQuadMesh(3, 3);
            auto quality = MeshQuality::fromMesh(mesh);

            EXPECT_EQ(quality.cellAspectRatio.size(), mesh.nCells);
        }

        TEST(MeshQualityTest, AspectRatioMinOne)
        {
            auto mesh = PolyMesh::createStructuredQuadMesh(3, 3);
            auto quality = MeshQuality::fromMesh(mesh);

            for (double ar : quality.cellAspectRatio)
            {
                EXPECT_GE(ar, 1.0);
            }
        }

        TEST(MeshQualityTest, UniformMeshUnitAspectRatio)
        {
            auto mesh = PolyMesh::createStructuredQuadMesh(5, 5);
            auto quality = MeshQuality::fromMesh(mesh);

            // A uniform square mesh should have aspect ratio close to 1
            for (double ar : quality.cellAspectRatio)
            {
                EXPECT_NEAR(ar, 1.0, 0.01);
            }
        }

        TEST(MeshQualityTest, NonOrthogonalitySize)
        {
            auto mesh = PolyMesh::createStructuredQuadMesh(3, 3);
            auto quality = MeshQuality::fromMesh(mesh);

            EXPECT_EQ(quality.cellNonOrthogonality.size(), mesh.nCells);
        }

        TEST(MeshQualityTest, NonOrthogonalityRange)
        {
            auto mesh = PolyMesh::createStructuredQuadMesh(3, 3);
            auto quality = MeshQuality::fromMesh(mesh);

            for (double nonOrth : quality.cellNonOrthogonality)
            {
                EXPECT_GE(nonOrth, 0.0);
                EXPECT_LE(nonOrth, 90.0);
            }
        }

        TEST(MeshQualityTest, UniformMeshVolumeRatioOne)
        {
            auto mesh = PolyMesh::createStructuredQuadMesh(5, 5);
            auto quality = MeshQuality::fromMesh(mesh);

            // All cells should have the same size
            EXPECT_NEAR(quality.minMaxVolumeRatio, 1.0, 1e-10);
        }

        // =============================================================================
        // Move Semantics Tests
        // =============================================================================

        TEST(PolyMeshMoveSemantics, MoveConstruct)
        {
            auto mesh1 = PolyMesh::createStructuredQuadMesh(3, 3);
            std::size_t originalNCells = mesh1.nCells;
            std::size_t originalNNodes = mesh1.nNodes;

            PolyMesh mesh2(std::move(mesh1));

            EXPECT_EQ(mesh2.nCells, originalNCells);
            EXPECT_EQ(mesh2.nNodes, originalNNodes);
            EXPECT_TRUE(mesh2.isAnalyzed());
        }

        TEST(PolyMeshMoveSemantics, MoveAssign)
        {
            auto mesh1 = PolyMesh::createStructuredQuadMesh(3, 3);
            std::size_t originalNCells = mesh1.nCells;

            PolyMesh mesh2;
            mesh2 = std::move(mesh1);

            EXPECT_EQ(mesh2.nCells, originalNCells);
            EXPECT_TRUE(mesh2.isAnalyzed());
        }

        // =============================================================================
        // Clear and Re-analyze Tests
        // =============================================================================

        TEST(PolyMeshAnalysis, ClearDerivedData)
        {
            auto mesh = PolyMesh::createStructuredQuadMesh(3, 3);
            EXPECT_TRUE(mesh.isAnalyzed());

            mesh.clearDerivedData();
            EXPECT_FALSE(mesh.isAnalyzed());
        }

        TEST(PolyMeshAnalysis, ReAnalyze)
        {
            auto mesh = PolyMesh::createStructuredQuadMesh(3, 3);
            auto originalCentroids = mesh.cellCentroids;

            mesh.clearDerivedData();
            mesh.analyzeMesh();

            EXPECT_TRUE(mesh.isAnalyzed());
            EXPECT_EQ(mesh.cellCentroids.size(), originalCentroids.size());

            // Centroids should be the same after re-analysis
            for (std::size_t i = 0; i < mesh.cellCentroids.size(); ++i)
            {
                EXPECT_NEAR(mesh.cellCentroids[i][0], originalCentroids[i][0], 1e-10);
                EXPECT_NEAR(mesh.cellCentroids[i][1], originalCentroids[i][1], 1e-10);
            }
        }

        // =============================================================================
        // Partition Tests (using hierarchical which doesn't require METIS)
        // =============================================================================

        TEST(PartitionTest, HierarchicalPartitionTwoParts)
        {
            auto mesh = PolyMesh::createStructuredQuadMesh(4, 4);
            auto parts = partitionWithHierarchical(mesh, 2);

            EXPECT_EQ(parts.size(), mesh.nCells);

            // Check that both partitions are used
            bool hasPart0 = false;
            bool hasPart1 = false;
            for (int p : parts)
            {
                if (p == 0)
                    hasPart0 = true;
                if (p == 1)
                    hasPart1 = true;
            }
            EXPECT_TRUE(hasPart0);
            EXPECT_TRUE(hasPart1);
        }

        TEST(PartitionTest, HierarchicalPartitionFourParts)
        {
            auto mesh = PolyMesh::createStructuredQuadMesh(4, 4);
            auto parts = partitionWithHierarchical(mesh, 4);

            EXPECT_EQ(parts.size(), mesh.nCells);

            // Check partition IDs are in valid range
            for (int p : parts)
            {
                EXPECT_GE(p, 0);
                EXPECT_LT(p, 4);
            }
        }

        TEST(PartitionTest, PartitionMeshFunctionHierarchical)
        {
            auto mesh = PolyMesh::createStructuredQuadMesh(4, 4);
            auto parts = partitionMesh(mesh, 2, "hierarchical");

            EXPECT_EQ(parts.size(), mesh.nCells);
        }

        TEST(PartitionTest, AdjacencyList)
        {
            auto mesh = PolyMesh::createStructuredQuadMesh(3, 3);
            auto adjList = getAdjacencyList(mesh);

            EXPECT_EQ(adjList.size(), mesh.nCells);

            // Interior cell should have 4 neighbors
            // Center cell in 3x3 grid is cell 4
            EXPECT_EQ(adjList[4].size(), 4u);

            // Corner cell should have 2 neighbors
            EXPECT_EQ(adjList[0].size(), 2u);
        }

    } // namespace testing
} // namespace fvm
