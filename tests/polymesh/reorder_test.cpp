#include <gtest/gtest.h>
#include "polymesh/poly_mesh.hpp"
#include "polymesh/reorder.hpp"
#include <algorithm>
#include <numeric>
#include <set>

namespace fvm
{
    namespace testing
    {

        // =============================================================================
        // Helper Functions
        // =============================================================================

        bool isValidPermutation(const std::vector<Index> &order, Index size)
        {
            if (order.size() != size)
                return false;

            std::vector<bool> seen(size, false);
            for (auto idx : order)
            {
                if (idx >= size || seen[idx])
                    return false;
                seen[idx] = true;
            }
            return true;
        }

        // =============================================================================
        // SparseMatrix Tests
        // =============================================================================

        TEST(SparseMatrixTest, DefaultConstruction)
        {
            SparseMatrix mat;
            EXPECT_EQ(mat.nRows, 0u);
        }

        TEST(SparseMatrixTest, ConstructionWithSize)
        {
            SparseMatrix mat(5);
            EXPECT_EQ(mat.nRows, 5u);
            EXPECT_EQ(mat.rowPtr.size(), 6u);
        }

        TEST(SparseMatrixTest, DegreeForEmptyRows)
        {
            SparseMatrix mat(3);
            // All rows are empty by default
            EXPECT_EQ(mat.degree(0), 0u);
            EXPECT_EQ(mat.degree(1), 0u);
            EXPECT_EQ(mat.degree(2), 0u);
        }

        // =============================================================================
        // Build Adjacency Tests
        // =============================================================================

        TEST(BuildAdjacencyTest, CellAdjacencySize)
        {
            auto mesh = PolyMesh::createStructuredQuadMesh(3, 3);
            auto adj = buildCellAdjacency(mesh);

            EXPECT_EQ(adj.nRows, mesh.nCells);
        }

        TEST(BuildAdjacencyTest, CellAdjacencyDegrees)
        {
            auto mesh = PolyMesh::createStructuredQuadMesh(3, 3);
            auto adj = buildCellAdjacency(mesh);

            // Corner cell should have degree 2
            EXPECT_EQ(adj.degree(0), 2u);

            // Interior cell should have degree 4
            EXPECT_EQ(adj.degree(4), 4u);
        }

        TEST(BuildAdjacencyTest, NodeAdjacencySize)
        {
            auto mesh = PolyMesh::createStructuredQuadMesh(3, 3);
            auto adj = buildNodeAdjacency(mesh);

            EXPECT_EQ(adj.nRows, mesh.nNodes);
        }

        TEST(BuildAdjacencyTest, CellAdjacencySubset)
        {
            auto mesh = PolyMesh::createStructuredQuadMesh(4, 4);
            auto adj = buildCellAdjacency(mesh, 8); // Only first 8 cells

            EXPECT_EQ(adj.nRows, 8u);
        }

        // =============================================================================
        // RCM Strategy Tests
        // =============================================================================

        TEST(RCMStrategyTest, ReturnsValidPermutation)
        {
            auto mesh = PolyMesh::createStructuredQuadMesh(4, 4);
            RCMStrategy strategy;

            auto order = strategy.getOrder(mesh);

            EXPECT_TRUE(isValidPermutation(order, mesh.nCells));
        }

        TEST(RCMStrategyTest, PermutationCoversAllCells)
        {
            auto mesh = PolyMesh::createStructuredQuadMesh(4, 4);
            RCMStrategy strategy;

            auto order = strategy.getOrder(mesh);

            std::set<std::size_t> uniqueIndices(order.begin(), order.end());
            EXPECT_EQ(uniqueIndices.size(), mesh.nCells);
        }

        // =============================================================================
        // GPS Strategy Tests
        // =============================================================================

        TEST(GPSStrategyTest, ReturnsValidPermutation)
        {
            auto mesh = PolyMesh::createStructuredQuadMesh(4, 4);
            GPSStrategy strategy;

            auto order = strategy.getOrder(mesh);

            EXPECT_TRUE(isValidPermutation(order, mesh.nCells));
        }

        // =============================================================================
        // Sloan Strategy Tests
        // =============================================================================

        TEST(SloanStrategyTest, ReturnsValidPermutation)
        {
            auto mesh = PolyMesh::createStructuredQuadMesh(4, 4);
            SloanStrategy strategy;

            auto order = strategy.getOrder(mesh);

            EXPECT_TRUE(isValidPermutation(order, mesh.nCells));
        }

        // =============================================================================
        // Spectral Strategy Tests
        // =============================================================================

        TEST(SpectralStrategyTest, ReturnsValidPermutation)
        {
            auto mesh = PolyMesh::createStructuredQuadMesh(4, 4);
            SpectralStrategy strategy;

            auto order = strategy.getOrder(mesh);

            EXPECT_TRUE(isValidPermutation(order, mesh.nCells));
        }

        // =============================================================================
        // Spatial Strategy Tests
        // =============================================================================

        TEST(SpatialXStrategyTest, ReturnsValidPermutation)
        {
            auto mesh = PolyMesh::createStructuredQuadMesh(4, 4);
            SpatialXStrategy strategy;

            auto order = strategy.getOrder(mesh);

            EXPECT_TRUE(isValidPermutation(order, mesh.nCells));
        }

        TEST(SpatialYStrategyTest, ReturnsValidPermutation)
        {
            auto mesh = PolyMesh::createStructuredQuadMesh(4, 4);
            SpatialYStrategy strategy;

            auto order = strategy.getOrder(mesh);

            EXPECT_TRUE(isValidPermutation(order, mesh.nCells));
        }

        // =============================================================================
        // Random Strategy Tests
        // =============================================================================

        TEST(RandomStrategyTest, ReturnsValidPermutation)
        {
            auto mesh = PolyMesh::createStructuredQuadMesh(4, 4);
            RandomStrategy strategy;

            auto order = strategy.getOrder(mesh);

            EXPECT_TRUE(isValidPermutation(order, mesh.nCells));
        }

        // =============================================================================
        // Create Strategy Factory Tests
        // =============================================================================

        TEST(CreateStrategyTest, CreateRCM)
        {
            auto strategy = createCellReorderStrategy("rcm");
            EXPECT_NE(strategy, nullptr);
        }

        TEST(CreateStrategyTest, CreateGPS)
        {
            auto strategy = createCellReorderStrategy("gps");
            EXPECT_NE(strategy, nullptr);
        }

        TEST(CreateStrategyTest, CreateSloan)
        {
            auto strategy = createCellReorderStrategy("sloan");
            EXPECT_NE(strategy, nullptr);
        }

        TEST(CreateStrategyTest, CreateSpectral)
        {
            auto strategy = createCellReorderStrategy("spectral");
            EXPECT_NE(strategy, nullptr);
        }

        TEST(CreateStrategyTest, CreateSpatialX)
        {
            auto strategy = createCellReorderStrategy("spatial_x");
            EXPECT_NE(strategy, nullptr);
        }

        TEST(CreateStrategyTest, CreateSpatialY)
        {
            auto strategy = createCellReorderStrategy("spatial_y");
            EXPECT_NE(strategy, nullptr);
        }

        TEST(CreateStrategyTest, CreateRandom)
        {
            auto strategy = createCellReorderStrategy("random");
            EXPECT_NE(strategy, nullptr);
        }

        TEST(CreateStrategyTest, UnknownStrategyReturnsNull)
        {
            EXPECT_THROW(createCellReorderStrategy("unknown_strategy"), std::invalid_argument);
        }

        // =============================================================================
        // RenumberCells Tests
        // =============================================================================

        TEST(RenumberCellsTest, RCMReordering)
        {
            auto mesh = PolyMesh::createStructuredQuadMesh(4, 4);
            std::size_t originalNCells = mesh.nCells;

            renumberCells(mesh, "rcm");

            EXPECT_EQ(mesh.nCells, originalNCells);
            EXPECT_EQ(mesh.cellNodeConnectivity.size(), originalNCells);
        }

        TEST(RenumberCellsTest, SpatialXReordering)
        {
            auto mesh = PolyMesh::createStructuredQuadMesh(4, 4);
            std::size_t originalNCells = mesh.nCells;

            renumberCells(mesh, "spatial_x");

            EXPECT_EQ(mesh.nCells, originalNCells);
        }

        TEST(RenumberCellsTest, SpatialYReordering)
        {
            auto mesh = PolyMesh::createStructuredQuadMesh(4, 4);
            std::size_t originalNCells = mesh.nCells;

            renumberCells(mesh, "spatial_y");

            EXPECT_EQ(mesh.nCells, originalNCells);
        }

        TEST(RenumberCellsTest, PartialReorder)
        {
            auto mesh = PolyMesh::createStructuredQuadMesh(4, 4);
            std::size_t originalNCells = mesh.nCells;

            // Reorder only the first 8 cells
            renumberCells(mesh, "rcm", 8);

            EXPECT_EQ(mesh.nCells, originalNCells);
        }

        TEST(RenumberCellsTest, InvalidStrategyNoChange)
        {
            auto mesh = PolyMesh::createStructuredQuadMesh(4, 4);
            auto originalConnectivity = mesh.cellNodeConnectivity;

            EXPECT_THROW(renumberCells(mesh, "invalid_strategy"), std::invalid_argument);

            // // Mesh should be unchanged with invalid strategy
            // EXPECT_EQ(mesh.cellNodeConnectivity.size(), originalConnectivity.size());
        }

        // =============================================================================
        // RenumberNodes Tests
        // =============================================================================

        TEST(RenumberNodesTest, RCMReordering)
        {
            auto mesh = PolyMesh::createStructuredQuadMesh(4, 4);
            std::size_t originalNNodes = mesh.nNodes;

            renumberNodes(mesh, "rcm");

            EXPECT_EQ(mesh.nNodes, originalNNodes);
            EXPECT_EQ(mesh.nodeCoords.size(), originalNNodes);
        }

        TEST(RenumberNodesTest, SequentialReordering)
        {
            auto mesh = PolyMesh::createStructuredQuadMesh(4, 4);
            std::size_t originalNNodes = mesh.nNodes;

            renumberNodes(mesh, "sequential");

            EXPECT_EQ(mesh.nNodes, originalNNodes);
        }

        TEST(RenumberNodesTest, ReverseReordering)
        {
            auto mesh = PolyMesh::createStructuredQuadMesh(4, 4);
            std::size_t originalNNodes = mesh.nNodes;

            renumberNodes(mesh, "reverse");

            EXPECT_EQ(mesh.nNodes, originalNNodes);
        }

        TEST(RenumberNodesTest, SpatialXReordering)
        {
            auto mesh = PolyMesh::createStructuredQuadMesh(4, 4);
            std::size_t originalNNodes = mesh.nNodes;

            renumberNodes(mesh, "spatial_x");

            EXPECT_EQ(mesh.nNodes, originalNNodes);
        }

        TEST(RenumberNodesTest, SpatialYReordering)
        {
            auto mesh = PolyMesh::createStructuredQuadMesh(4, 4);
            std::size_t originalNNodes = mesh.nNodes;

            renumberNodes(mesh, "spatial_y");

            EXPECT_EQ(mesh.nNodes, originalNNodes);
        }

        TEST(RenumberNodesTest, RandomReordering)
        {
            auto mesh = PolyMesh::createStructuredQuadMesh(4, 4);
            std::size_t originalNNodes = mesh.nNodes;

            renumberNodes(mesh, "random");

            EXPECT_EQ(mesh.nNodes, originalNNodes);
        }

        // =============================================================================
        // Mesh Integrity After Reordering Tests
        // =============================================================================

        TEST(ReorderIntegrityTest, CellVolumesPreserved)
        {
            auto mesh = PolyMesh::createStructuredQuadMesh(4, 4);

            double totalVolumeBefore = 0.0;
            for (double vol : mesh.cellVolumes)
            {
                totalVolumeBefore += vol;
            }

            renumberCells(mesh, "rcm");
            mesh.analyzeMesh();

            double totalVolumeAfter = 0.0;
            for (double vol : mesh.cellVolumes)
            {
                totalVolumeAfter += vol;
            }

            EXPECT_NEAR(totalVolumeBefore, totalVolumeAfter, 1e-10);
        }

        TEST(ReorderIntegrityTest, NodeCountPreserved)
        {
            auto mesh = PolyMesh::createStructuredQuadMesh(4, 4);
            std::size_t originalNNodes = mesh.nNodes;

            renumberNodes(mesh, "rcm");

            EXPECT_EQ(mesh.nNodes, originalNNodes);
        }

        TEST(ReorderIntegrityTest, CellCountPreserved)
        {
            auto mesh = PolyMesh::createStructuredQuadMesh(4, 4);
            std::size_t originalNCells = mesh.nCells;

            renumberCells(mesh, "spectral");

            EXPECT_EQ(mesh.nCells, originalNCells);
        }

        TEST(ReorderIntegrityTest, ConnectivityIndicesValid)
        {
            auto mesh = PolyMesh::createStructuredQuadMesh(4, 4);

            renumberCells(mesh, "rcm");
            renumberNodes(mesh, "rcm");

            for (const auto &cell : mesh.cellNodeConnectivity)
            {
                for (auto nodeIdx : cell)
                {
                    EXPECT_LT(nodeIdx, mesh.nNodes);
                }
            }
        }

    } // namespace testing
} // namespace fvm
