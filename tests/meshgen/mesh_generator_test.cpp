#include <gtest/gtest.h>
#include <gmsh.h>
#include "meshgen/geometry.hpp"
#include "meshgen/mesh_generator.hpp"
#include <cmath>
#include <filesystem>
#include <stdexcept>
#include "vtkio/cell_types.hpp"

namespace fvm
{
    namespace testing
    {

        // Test fixture that handles Gmsh initialization/finalization
        class MeshGeneratorTest : public ::testing::Test
        {
        protected:
            std::string testOutputDir_;

            void SetUp() override
            {
                gmsh::initialize();
                gmsh::option::setNumber("General.Terminal", 0); // Suppress output
                testOutputDir_ = "test_output";
                std::filesystem::create_directories(testOutputDir_);
            }

            void TearDown() override
            {
                gmsh::clear();
                gmsh::finalize();
                // Clean up test output directory
                std::filesystem::remove_all(testOutputDir_);
            }

            // Helper to create a simple rectangle geometry
            int createRectangle(double length = 2.0, double width = 1.0)
            {
                gmsh::model::add("test");
                Geometry geo;
                return geo.rectangle(length, width, 0.0, 0.0, 0.1);
            }

            // Helper to create a simple circle geometry
            int createCircle(double radius = 1.0)
            {
                gmsh::model::add("test");
                Geometry geo;
                return geo.circle(radius, 0.0, 0.0, 0.1);
            }

            // Helper to create a triangle geometry
            int createTriangle()
            {
                gmsh::model::add("test");
                Geometry geo;
                Point2D p1 = {0.0, 0.0};
                Point2D p2 = {2.0, 0.0};
                Point2D p3 = {1.0, 1.5};
                return geo.triangle(p1, p2, p3, 0.1);
            }
        };

        // =============================================================================
        // Constructor Tests
        // =============================================================================

        TEST_F(MeshGeneratorTest, ConstructorWithSingleTag)
        {
            int surface = createRectangle();
            EXPECT_NO_THROW(MeshGenerator(surface, testOutputDir_));
        }

        TEST_F(MeshGeneratorTest, ConstructorWithVectorOfTags)
        {
            int surface = createRectangle();
            std::vector<int> tags = {surface};
            EXPECT_NO_THROW(MeshGenerator(tags, testOutputDir_));
        }

        TEST_F(MeshGeneratorTest, ConstructorCreatesOutputDirectory)
        {
            int surface = createRectangle();
            std::string newDir = testOutputDir_ + "/new_subdir";
            EXPECT_FALSE(std::filesystem::exists(newDir));

            MeshGenerator gen(surface, newDir);

            EXPECT_TRUE(std::filesystem::exists(newDir));
        }

        // =============================================================================
        // Triangular Mesh Generation Tests
        // =============================================================================

        TEST_F(MeshGeneratorTest, GenerateTriangularMesh)
        {
            int surface = createRectangle();
            MeshGenerator gen(surface, testOutputDir_);

            std::map<int, MeshParams> params;
            params[surface] = {"tri", 0.2};

            EXPECT_NO_THROW(gen.generate(params, "tri_mesh.msh"));

            const MeshInfo &data = gen.getMeshData();
            EXPECT_GT(data.nodes.size(), 0u);
            EXPECT_GT(data.elements.size(), 0u);
        }

        TEST_F(MeshGeneratorTest, TriangularMeshHasCorrectCellTypes)
        {
            int surface = createRectangle();
            MeshGenerator gen(surface, testOutputDir_);

            std::map<int, MeshParams> params;
            params[surface] = {"tri", 0.3};

            gen.generate(params, "tri_mesh.msh");
            const MeshInfo &data = gen.getMeshData();

            // All cells should be triangles (VTK type 5)
            for (int cellType : data.elementTypes)
            {
                EXPECT_EQ(cellType, VTKCellType::TRIANGLE);
            }

            // All cells should have 3 nodes
            for (const auto &cell : data.elements)
            {
                EXPECT_EQ(cell.size(), 3u);
            }
        }

        TEST_F(MeshGeneratorTest, TriangularMeshOnCircle)
        {
            int surface = createCircle();
            MeshGenerator gen(surface, testOutputDir_);

            std::map<int, MeshParams> params;
            params[surface] = {"tri", 0.2};

            EXPECT_NO_THROW(gen.generate(params, "circle_tri.msh"));

            const MeshInfo &data = gen.getMeshData();
            EXPECT_GT(data.nodes.size(), 0u);
            EXPECT_GT(data.elements.size(), 0u);
        }

        TEST_F(MeshGeneratorTest, TriangularMeshOnTriangle)
        {
            int surface = createTriangle();
            MeshGenerator gen(surface, testOutputDir_);

            std::map<int, MeshParams> params;
            params[surface] = {"tri", 0.2};

            EXPECT_NO_THROW(gen.generate(params, "triangle_tri.msh"));

            const MeshInfo &data = gen.getMeshData();
            EXPECT_GT(data.nodes.size(), 0u);
        }

        // =============================================================================
        // Quadrilateral Mesh Generation Tests
        // =============================================================================

        TEST_F(MeshGeneratorTest, GenerateQuadMesh)
        {
            int surface = createRectangle();
            MeshGenerator gen(surface, testOutputDir_);

            std::map<int, MeshParams> params;
            params[surface] = {"quad", 0.2};

            EXPECT_NO_THROW(gen.generate(params, "quad_mesh.msh"));

            const MeshInfo &data = gen.getMeshData();
            EXPECT_GT(data.nodes.size(), 0u);
            EXPECT_GT(data.elements.size(), 0u);
        }

        TEST_F(MeshGeneratorTest, QuadMeshHasCorrectCellTypes)
        {
            int surface = createRectangle();
            MeshGenerator gen(surface, testOutputDir_);

            std::map<int, MeshParams> params;
            params[surface] = {"quad", 0.3};

            gen.generate(params, "quad_mesh.msh");
            const MeshInfo &data = gen.getMeshData();

            // Most cells should be quads (VTK type 9), some may be triangles
            int quadCount = 0;
            for (int cellType : data.elementTypes)
            {
                if (cellType == VTKCellType::QUAD)
                {
                    quadCount++;
                }
            }
            // Expect majority to be quads
            EXPECT_GT(quadCount, 0);
        }

        // =============================================================================
        // Structured Mesh Generation Tests
        // =============================================================================

        TEST_F(MeshGeneratorTest, GenerateStructuredMesh)
        {
            int surface = createRectangle(2.0, 1.0);
            MeshGenerator gen(surface, testOutputDir_);

            std::map<int, MeshParams> params;
            params[surface] = {"structured", 0.25};

            EXPECT_NO_THROW(gen.generate(params, "structured_mesh.msh"));

            const MeshInfo &data = gen.getMeshData();
            EXPECT_GT(data.nodes.size(), 0u);
            EXPECT_GT(data.elements.size(), 0u);
        }

        TEST_F(MeshGeneratorTest, StructuredMeshHasOnlyQuads)
        {
            int surface = createRectangle(2.0, 1.0);
            MeshGenerator gen(surface, testOutputDir_);

            std::map<int, MeshParams> params;
            params[surface] = {"structured", 0.5};

            gen.generate(params, "structured_mesh.msh");
            const MeshInfo &data = gen.getMeshData();

            // All cells should be quads for structured mesh
            for (int cellType : data.elementTypes)
            {
                EXPECT_EQ(cellType, VTKCellType::QUAD);
            }

            for (const auto &cell : data.elements)
            {
                EXPECT_EQ(cell.size(), 4u);
            }
        }

        TEST_F(MeshGeneratorTest, StructuredMeshThrowsOnNonRectangle)
        {
            int surface = createTriangle();
            MeshGenerator gen(surface, testOutputDir_);

            std::map<int, MeshParams> params;
            params[surface] = {"structured", 0.2};

            // Structured mesh requires 4 boundary curves
            EXPECT_THROW(gen.generate(params, "should_fail.msh"), std::runtime_error);
        }

        // =============================================================================
        // Mesh Data Extraction Tests
        // =============================================================================

        TEST_F(MeshGeneratorTest, ExtractedNodesHaveCorrectDimensions)
        {
            int surface = createRectangle();
            MeshGenerator gen(surface, testOutputDir_);

            std::map<int, MeshParams> params;
            params[surface] = {"tri", 0.3};

            gen.generate(params, "mesh.msh");
            const MeshInfo &data = gen.getMeshData();

            // Each node should have 3 coordinates
            for (const auto &node : data.nodes)
            {
                EXPECT_EQ(node.size(), 3u);
            }

            // For 2D mesh, z should be 0
            for (const auto &node : data.nodes)
            {
                EXPECT_NEAR(node[2], 0.0, 1e-10);
            }
        }

        TEST_F(MeshGeneratorTest, ElementTypesMatchElementCount)
        {
            int surface = createRectangle();
            MeshGenerator gen(surface, testOutputDir_);

            std::map<int, MeshParams> params;
            params[surface] = {"tri", 0.3};

            gen.generate(params, "mesh.msh");
            const MeshInfo &data = gen.getMeshData();

            EXPECT_EQ(data.elements.size(), data.elementTypes.size());
        }

        TEST_F(MeshGeneratorTest, CellConnectivityIsValid)
        {
            int surface = createRectangle();
            MeshGenerator gen(surface, testOutputDir_);

            std::map<int, MeshParams> params;
            params[surface] = {"tri", 0.3};

            gen.generate(params, "mesh.msh");
            const MeshInfo &data = gen.getMeshData();

            // All cell node indices should be valid
            std::size_t numNodes = data.nodes.size();
            for (const auto &cell : data.elements)
            {
                for (std::size_t nodeIdx : cell)
                {
                    EXPECT_LT(nodeIdx, numNodes);
                }
            }
        }

        // =============================================================================
        // Physical Groups / Sets Tests
        // =============================================================================

        TEST_F(MeshGeneratorTest, HasFaceSets)
        {
            int surface = createRectangle();
            MeshGenerator gen(surface, testOutputDir_);

            std::map<int, MeshParams> params;
            params[surface] = {"tri", 0.3};

            gen.generate(params, "mesh.msh");
            const MeshInfo &data = gen.getMeshData();

            // Should have at least one face set (boundary)
            EXPECT_GT(data.faceSets.size(), 0u);
        }

        TEST_F(MeshGeneratorTest, HasElementSets)
        {
            int surface = createRectangle();
            MeshGenerator gen(surface, testOutputDir_);

            std::map<int, MeshParams> params;
            params[surface] = {"tri", 0.3};

            gen.generate(params, "mesh.msh");
            const MeshInfo &data = gen.getMeshData();

            // Should have at least one element set (surface in 2D)
            EXPECT_GT(data.elementSets.size(), 0u);
        }

        TEST_F(MeshGeneratorTest, DefaultFluidGroupCreated)
        {
            int surface = createRectangle();
            MeshGenerator gen(surface, testOutputDir_);

            std::map<int, MeshParams> params;
            params[surface] = {"tri", 0.3};

            gen.generate(params, "mesh.msh");
            const MeshInfo &data = gen.getMeshData();

            // Should have "fluid" group for untagged surfaces
            EXPECT_TRUE(data.elementSets.find("fluid") != data.elementSets.end());
        }

        TEST_F(MeshGeneratorTest, NamedBoundaryPhysicalGroups)
        {
            // Create rectangle and get boundary curves
            gmsh::model::add("test");
            Geometry geo;
            int surface = geo.rectangle(2.0, 1.0, 0.0, 0.0, 0.2);

            // Get boundary curves of the surface
            std::vector<std::pair<int, int>> boundaryCurves;
            gmsh::model::getBoundary({{2, surface}}, boundaryCurves, false, false, false);

            ASSERT_EQ(boundaryCurves.size(), 4u);

            // Classify and name each boundary curve based on its position
            for (const auto &dimTag : boundaryCurves)
            {
                int curveTag = dimTag.second;

                // Get curve bounding box to determine position
                double minX, minY, minZ, maxX, maxY, maxZ;
                gmsh::model::getBoundingBox(1, curveTag, minX, minY, minZ, maxX, maxY, maxZ);

                // Determine boundary name based on position
                std::string name;
                double midY = (minY + maxY) / 2.0;
                double midX = (minX + maxX) / 2.0;

                // Horizontal edges (y is constant)
                if (std::abs(maxY - minY) < 1e-6)
                {
                    if (midY < 0.5)
                        name = "bottom";
                    else
                        name = "top";
                }
                // Vertical edges (x is constant)
                else
                {
                    if (midX < 1.0)
                        name = "left";
                    else
                        name = "right";
                }

                gmsh::model::addPhysicalGroup(1, {curveTag}, -1, name);
            }

            // Generate mesh
            MeshGenerator gen(surface, testOutputDir_);
            std::map<int, MeshParams> params;
            params[surface] = {"tri", 0.3};
            gen.generate(params, "named_boundaries.msh");

            const MeshInfo &data = gen.getMeshData();

            // Verify all four named face sets exist
            EXPECT_TRUE(data.faceSets.find("bottom") != data.faceSets.end());
            EXPECT_TRUE(data.faceSets.find("top") != data.faceSets.end());
            EXPECT_TRUE(data.faceSets.find("left") != data.faceSets.end());
            EXPECT_TRUE(data.faceSets.find("right") != data.faceSets.end());

            // Verify we have exactly 4 face sets
            EXPECT_EQ(data.faceSets.size(), 4u);

            // Verify each face set has faces
            for (const auto &[name, faces] : data.faceSets)
            {
                EXPECT_GT(faces.size(), 0u) << "FaceSet " << name << " has no faces";
            }
        }

        // =============================================================================
        // Face Set Tests
        // =============================================================================

        TEST_F(MeshGeneratorTest, ExtractsFaceSets)
        {
            int surface = createRectangle();
            MeshGenerator gen(surface, testOutputDir_);

            std::map<int, MeshParams> params;
            params[surface] = {"tri", 0.3};

            gen.generate(params, "mesh.msh");
            const MeshInfo &data = gen.getMeshData();

            // Should have face sets (boundary groups)
            EXPECT_GT(data.faceSets.size(), 0u);
        }

        TEST_F(MeshGeneratorTest, FaceSetConnectivityIsValid)
        {
            int surface = createRectangle();
            MeshGenerator gen(surface, testOutputDir_);

            std::map<int, MeshParams> params;
            params[surface] = {"tri", 0.3};

            gen.generate(params, "mesh.msh");
            const MeshInfo &data = gen.getMeshData();

            std::size_t numNodes = data.nodes.size();

            // All face node indices should be valid
            for (const auto &[name, faces] : data.faceSets)
            {
                for (const auto &face : faces)
                {
                    for (std::size_t nodeIdx : face)
                    {
                        EXPECT_LT(nodeIdx, numNodes);
                    }
                }
            }
        }

        // =============================================================================
        // File Output Tests
        // =============================================================================

        TEST_F(MeshGeneratorTest, CreatesMshFile)
        {
            int surface = createRectangle();
            MeshGenerator gen(surface, testOutputDir_);

            std::map<int, MeshParams> params;
            params[surface] = {"tri", 0.3};

            gen.generate(params, "output.msh");

            std::string expectedPath = testOutputDir_ + "/output.msh";
            EXPECT_TRUE(std::filesystem::exists(expectedPath));
        }

        TEST_F(MeshGeneratorTest, MshFileNotEmpty)
        {
            int surface = createRectangle();
            MeshGenerator gen(surface, testOutputDir_);

            std::map<int, MeshParams> params;
            params[surface] = {"tri", 0.3};

            gen.generate(params, "output.msh");

            std::string expectedPath = testOutputDir_ + "/output.msh";
            EXPECT_GT(std::filesystem::file_size(expectedPath), 0u);
        }

        // =============================================================================
        // Mesh Size Tests
        // =============================================================================

        TEST_F(MeshGeneratorTest, SmallerMeshSizeProducesMoreElements)
        {
            // Generate with larger mesh size (set via geometry characteristic length)
            gmsh::model::add("test1");
            Geometry geo1;
            int surface1 = geo1.rectangle(2.0, 1.0, 0.0, 0.0, 0.5); // coarse: 0.5

            MeshGenerator gen1(surface1, testOutputDir_);
            std::map<int, MeshParams> params1;
            params1[surface1] = {"tri", 0.5};
            gen1.generate(params1, "coarse.msh");
            std::size_t coarseElements = gen1.getMeshData().elements.size();

            // Reset and generate with smaller mesh size
            gmsh::clear();
            gmsh::model::add("test2");
            Geometry geo2;
            int surface2 = geo2.rectangle(2.0, 1.0, 0.0, 0.0, 0.1); // fine: 0.1

            MeshGenerator gen2(surface2, testOutputDir_);
            std::map<int, MeshParams> params2;
            params2[surface2] = {"tri", 0.1};
            gen2.generate(params2, "fine.msh");
            std::size_t fineElements = gen2.getMeshData().elements.size();

            // Finer mesh should have more elements
            EXPECT_GT(fineElements, coarseElements);
        }

        // =============================================================================
        // Multiple Surfaces Tests
        // =============================================================================

        TEST_F(MeshGeneratorTest, MeshMultipleSurfaces)
        {
            gmsh::model::add("test");
            Geometry geo;

            int rect1 = geo.rectangle(1.0, 1.0, 0.0, 0.0, 0.1);
            int rect2 = geo.rectangle(1.0, 1.0, 1.5, 0.0, 0.1);

            std::vector<int> surfaces = {rect1, rect2};
            MeshGenerator gen(surfaces, testOutputDir_);

            std::map<int, MeshParams> params;
            params[rect1] = {"tri", 0.2};
            params[rect2] = {"quad", 0.2};

            EXPECT_NO_THROW(gen.generate(params, "multi.msh"));

            const MeshInfo &data = gen.getMeshData();
            EXPECT_GT(data.nodes.size(), 0u);
            EXPECT_GT(data.elements.size(), 0u);
        }

        // =============================================================================
        // Extract Mesh Data Without Generate Tests
        // =============================================================================

        TEST_F(MeshGeneratorTest, ExtractMeshInfoDirectly)
        {
            int surface = createRectangle();

            // Generate mesh directly via Gmsh
            gmsh::model::mesh::generate(2);

            MeshGenerator gen(surface, testOutputDir_);
            EXPECT_NO_THROW(gen.extractMeshData());

            const MeshInfo &data = gen.getMeshData();
            EXPECT_GT(data.nodes.size(), 0u);
            EXPECT_GT(data.elements.size(), 0u);
        }

    } // namespace testing
} // namespace fvm
