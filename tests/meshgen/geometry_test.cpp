#include <gtest/gtest.h>
#include <gmsh.h>
#include "meshgen/geometry.hpp"
#include <cmath>
#include <stdexcept>

namespace fvm
{
    namespace testing
    {

        // Test fixture that handles Gmsh initialization/finalization
        class GeometryTest : public ::testing::Test
        {
        protected:
            void SetUp() override
            {
                gmsh::initialize();
                gmsh::option::setNumber("General.Terminal", 0); // Suppress output
            }

            void TearDown() override
            {
                gmsh::clear();
                gmsh::finalize();
            }
        };

        // =============================================================================
        // Constructor and getName Tests
        // =============================================================================

        TEST_F(GeometryTest, ConstructorWithDefaultName)
        {
            Geometry geo;
            EXPECT_EQ(geo.getName(), "Default Geometry");
        }

        TEST_F(GeometryTest, ConstructorWithEmptyName)
        {
            Geometry geo("");
            EXPECT_EQ(geo.getName(), "Default Geometry");
        }

        TEST_F(GeometryTest, ConstructorWithCustomName)
        {
            Geometry geo("My Custom Geometry");
            EXPECT_EQ(geo.getName(), "My Custom Geometry");
        }

        // =============================================================================
        // Rectangle Tests
        // =============================================================================

        TEST_F(GeometryTest, RectangleCreatesValidSurface)
        {
            gmsh::model::add("test");
            Geometry geo("RectTest");

            int surfaceTag = geo.rectangle(2.0, 1.0, 0.0, 0.0, 0.1);

            EXPECT_GT(surfaceTag, 0);
        }

        TEST_F(GeometryTest, RectangleWithOffset)
        {
            gmsh::model::add("test");
            Geometry geo("RectOffsetTest");

            int surfaceTag = geo.rectangle(2.0, 1.0, 1.0, 2.0, 0.1);

            EXPECT_GT(surfaceTag, 0);

            BoundingBox bbox = geo.getBoundingBox();
            // Check bounding box matches expected dimensions
            EXPECT_NEAR(bbox[0], 1.0, 1e-10); // min_x
            EXPECT_NEAR(bbox[1], 2.0, 1e-10); // min_y
            EXPECT_NEAR(bbox[3], 3.0, 1e-10); // max_x = x + length
            EXPECT_NEAR(bbox[4], 3.0, 1e-10); // max_y = y + width
        }

        TEST_F(GeometryTest, RectangleBoundingBox)
        {
            gmsh::model::add("test");
            Geometry geo;

            geo.rectangle(4.0, 3.0, 0.0, 0.0, 0.1);
            BoundingBox bbox = geo.getBoundingBox();

            EXPECT_NEAR(bbox[0], 0.0, 1e-10); // min_x
            EXPECT_NEAR(bbox[1], 0.0, 1e-10); // min_y
            EXPECT_NEAR(bbox[2], 0.0, 1e-10); // min_z
            EXPECT_NEAR(bbox[3], 4.0, 1e-10); // max_x
            EXPECT_NEAR(bbox[4], 3.0, 1e-10); // max_y
            EXPECT_NEAR(bbox[5], 0.0, 1e-10); // max_z
        }

        // =============================================================================
        // Circle Tests
        // =============================================================================

        TEST_F(GeometryTest, CircleCreatesValidSurface)
        {
            gmsh::model::add("test");
            Geometry geo("CircleTest");

            int surfaceTag = geo.circle(1.0, 0.0, 0.0, 0.1);

            EXPECT_GT(surfaceTag, 0);
        }

        TEST_F(GeometryTest, CircleBoundingBox)
        {
            gmsh::model::add("test");
            Geometry geo;

            double radius = 2.0;
            double cx = 1.0, cy = 1.0;
            geo.circle(radius, cx, cy, 0.1);

            BoundingBox bbox = geo.getBoundingBox();

            EXPECT_NEAR(bbox[0], cx - radius, 1e-10); // min_x
            EXPECT_NEAR(bbox[1], cy - radius, 1e-10); // min_y
            EXPECT_NEAR(bbox[3], cx + radius, 1e-10); // max_x
            EXPECT_NEAR(bbox[4], cy + radius, 1e-10); // max_y
        }

        TEST_F(GeometryTest, CircleWithOffset)
        {
            gmsh::model::add("test");
            Geometry geo;

            double radius = 1.5;
            double cx = 5.0, cy = 3.0;
            int surfaceTag = geo.circle(radius, cx, cy, 0.1);

            EXPECT_GT(surfaceTag, 0);

            BoundingBox bbox = geo.getBoundingBox();
            EXPECT_NEAR(bbox[0], cx - radius, 1e-10);
            EXPECT_NEAR(bbox[1], cy - radius, 1e-10);
            EXPECT_NEAR(bbox[3], cx + radius, 1e-10);
            EXPECT_NEAR(bbox[4], cy + radius, 1e-10);
        }

        // =============================================================================
        // Triangle Tests
        // =============================================================================

        TEST_F(GeometryTest, TriangleCreatesValidSurface)
        {
            gmsh::model::add("test");
            Geometry geo("TriangleTest");

            Point2D p1 = {0.0, 0.0};
            Point2D p2 = {1.0, 0.0};
            Point2D p3 = {0.5, 1.0};

            int surfaceTag = geo.triangle(p1, p2, p3, 0.1);

            EXPECT_GT(surfaceTag, 0);
        }

        TEST_F(GeometryTest, TriangleBoundingBox)
        {
            gmsh::model::add("test");
            Geometry geo;

            Point2D p1 = {0.0, 0.0};
            Point2D p2 = {3.0, 0.0};
            Point2D p3 = {1.5, 2.0};

            geo.triangle(p1, p2, p3, 0.1);
            BoundingBox bbox = geo.getBoundingBox();

            EXPECT_NEAR(bbox[0], 0.0, 1e-10); // min_x
            EXPECT_NEAR(bbox[1], 0.0, 1e-10); // min_y
            EXPECT_NEAR(bbox[3], 3.0, 1e-10); // max_x
            EXPECT_NEAR(bbox[4], 2.0, 1e-10); // max_y
        }

        TEST_F(GeometryTest, RightTriangle)
        {
            gmsh::model::add("test");
            Geometry geo;

            Point2D p1 = {0.0, 0.0};
            Point2D p2 = {4.0, 0.0};
            Point2D p3 = {0.0, 3.0};

            int surfaceTag = geo.triangle(p1, p2, p3, 0.1);

            EXPECT_GT(surfaceTag, 0);

            BoundingBox bbox = geo.getBoundingBox();
            EXPECT_NEAR(bbox[0], 0.0, 1e-10);
            EXPECT_NEAR(bbox[1], 0.0, 1e-10);
            EXPECT_NEAR(bbox[3], 4.0, 1e-10);
            EXPECT_NEAR(bbox[4], 3.0, 1e-10);
        }

        // =============================================================================
        // Ellipse Tests
        // =============================================================================

        TEST_F(GeometryTest, EllipseCreatesValidSurface)
        {
            gmsh::model::add("test");
            Geometry geo("EllipseTest");

            int surfaceTag = geo.ellipse(2.0, 1.0, 0.0, 0.0, 0.1);

            EXPECT_GT(surfaceTag, 0);
        }

        TEST_F(GeometryTest, EllipseBoundingBox)
        {
            gmsh::model::add("test");
            Geometry geo;

            double r1 = 3.0, r2 = 2.0;
            double cx = 0.0, cy = 0.0;

            geo.ellipse(r1, r2, cx, cy, 0.1);
            BoundingBox bbox = geo.getBoundingBox();

            EXPECT_NEAR(bbox[0], -r1, 1e-10); // min_x
            EXPECT_NEAR(bbox[1], -r2, 1e-10); // min_y
            EXPECT_NEAR(bbox[3], r1, 1e-10);  // max_x
            EXPECT_NEAR(bbox[4], r2, 1e-10);  // max_y
        }

        TEST_F(GeometryTest, EllipseWithOffset)
        {
            gmsh::model::add("test");
            Geometry geo;

            double r1 = 2.0, r2 = 1.5;
            double cx = 5.0, cy = 3.0;

            geo.ellipse(r1, r2, cx, cy, 0.1);
            BoundingBox bbox = geo.getBoundingBox();

            EXPECT_NEAR(bbox[0], cx - r1, 1e-10);
            EXPECT_NEAR(bbox[1], cy - r2, 1e-10);
            EXPECT_NEAR(bbox[3], cx + r1, 1e-10);
            EXPECT_NEAR(bbox[4], cy + r2, 1e-10);
        }

        TEST_F(GeometryTest, CircularEllipse)
        {
            // When r1 == r2, ellipse should behave like a circle
            gmsh::model::add("test");
            Geometry geo;

            double r = 2.0;
            geo.ellipse(r, r, 0.0, 0.0, 0.1);
            BoundingBox bbox = geo.getBoundingBox();

            EXPECT_NEAR(bbox[0], -r, 1e-10);
            EXPECT_NEAR(bbox[1], -r, 1e-10);
            EXPECT_NEAR(bbox[3], r, 1e-10);
            EXPECT_NEAR(bbox[4], r, 1e-10);
        }

        // =============================================================================
        // Polygon Tests
        // =============================================================================

        TEST_F(GeometryTest, PolygonCreatesValidSurface)
        {
            gmsh::model::add("test");
            Geometry geo("PolygonTest");

            std::vector<Point2D> points = {
                {0.0, 0.0},
                {2.0, 0.0},
                {2.0, 1.0},
                {0.0, 1.0}};

            int surfaceTag = geo.polygon(points, false, 0.1);

            EXPECT_GT(surfaceTag, 0);
        }

        TEST_F(GeometryTest, PolygonWithThreePoints)
        {
            gmsh::model::add("test");
            Geometry geo;

            std::vector<Point2D> points = {
                {0.0, 0.0},
                {1.0, 0.0},
                {0.5, 1.0}};

            int surfaceTag = geo.polygon(points, false, 0.1);

            EXPECT_GT(surfaceTag, 0);
        }

        TEST_F(GeometryTest, PolygonBoundingBox)
        {
            gmsh::model::add("test");
            Geometry geo;

            std::vector<Point2D> points = {
                {1.0, 1.0},
                {4.0, 1.0},
                {4.0, 3.0},
                {2.5, 4.0},
                {1.0, 3.0}};

            geo.polygon(points, false, 0.1);
            BoundingBox bbox = geo.getBoundingBox();

            EXPECT_NEAR(bbox[0], 1.0, 1e-10); // min_x
            EXPECT_NEAR(bbox[1], 1.0, 1e-10); // min_y
            EXPECT_NEAR(bbox[3], 4.0, 1e-10); // max_x
            EXPECT_NEAR(bbox[4], 4.0, 1e-10); // max_y
        }

        TEST_F(GeometryTest, PolygonThrowsWithLessThanThreePoints)
        {
            gmsh::model::add("test");
            Geometry geo;

            std::vector<Point2D> twoPoints = {
                {0.0, 0.0},
                {1.0, 0.0}};

            EXPECT_THROW(geo.polygon(twoPoints, false, 0.1), std::invalid_argument);
        }

        TEST_F(GeometryTest, PolygonThrowsWithOnePoint)
        {
            gmsh::model::add("test");
            Geometry geo;

            std::vector<Point2D> onePoint = {{0.0, 0.0}};

            EXPECT_THROW(geo.polygon(onePoint, false, 0.1), std::invalid_argument);
        }

        TEST_F(GeometryTest, PolygonThrowsWithEmptyPoints)
        {
            gmsh::model::add("test");
            Geometry geo;

            std::vector<Point2D> emptyPoints;

            EXPECT_THROW(geo.polygon(emptyPoints, false, 0.1), std::invalid_argument);
        }

        // =============================================================================
        // Convex Hull Tests (via polygon method)
        // =============================================================================

        TEST_F(GeometryTest, ConvexHullSquarePoints)
        {
            gmsh::model::add("test");
            Geometry geo;

            // Points forming a square - convex hull should include all
            std::vector<Point2D> points = {
                {0.0, 0.0},
                {1.0, 0.0},
                {1.0, 1.0},
                {0.0, 1.0}};

            int surfaceTag = geo.polygon(points, true, 0.1);

            EXPECT_GT(surfaceTag, 0);

            BoundingBox bbox = geo.getBoundingBox();
            EXPECT_NEAR(bbox[0], 0.0, 1e-10);
            EXPECT_NEAR(bbox[1], 0.0, 1e-10);
            EXPECT_NEAR(bbox[3], 1.0, 1e-10);
            EXPECT_NEAR(bbox[4], 1.0, 1e-10);
        }

        TEST_F(GeometryTest, ConvexHullWithInteriorPoints)
        {
            gmsh::model::add("test");
            Geometry geo;

            // Square with interior point - convex hull should exclude interior point
            std::vector<Point2D> points = {
                {0.0, 0.0},
                {2.0, 0.0},
                {2.0, 2.0},
                {0.0, 2.0},
                {1.0, 1.0} // interior point
            };

            int surfaceTag = geo.polygon(points, true, 0.1);

            EXPECT_GT(surfaceTag, 0);

            // Bounding box should match the outer square
            BoundingBox bbox = geo.getBoundingBox();
            EXPECT_NEAR(bbox[0], 0.0, 1e-10);
            EXPECT_NEAR(bbox[1], 0.0, 1e-10);
            EXPECT_NEAR(bbox[3], 2.0, 1e-10);
            EXPECT_NEAR(bbox[4], 2.0, 1e-10);
        }

        TEST_F(GeometryTest, ConvexHullTriangle)
        {
            gmsh::model::add("test");
            Geometry geo;

            std::vector<Point2D> points = {
                {0.0, 0.0},
                {3.0, 0.0},
                {1.5, 2.5}};

            int surfaceTag = geo.polygon(points, true, 0.1);

            EXPECT_GT(surfaceTag, 0);
        }

        TEST_F(GeometryTest, ConvexHullRandomPoints)
        {
            gmsh::model::add("test");
            Geometry geo;

            // Points where convex hull should be the outer boundary
            std::vector<Point2D> points = {
                {0.0, 0.0},
                {4.0, 0.0},
                {4.0, 3.0},
                {0.0, 3.0},
                {1.0, 1.0}, // interior
                {2.0, 2.0}, // interior
                {3.0, 1.5}  // interior
            };

            int surfaceTag = geo.polygon(points, true, 0.1);

            EXPECT_GT(surfaceTag, 0);

            BoundingBox bbox = geo.getBoundingBox();
            EXPECT_NEAR(bbox[0], 0.0, 1e-10);
            EXPECT_NEAR(bbox[1], 0.0, 1e-10);
            EXPECT_NEAR(bbox[3], 4.0, 1e-10);
            EXPECT_NEAR(bbox[4], 3.0, 1e-10);
        }

        // =============================================================================
        // RectangleWithPartitions Tests
        // =============================================================================

        TEST_F(GeometryTest, RectangleWithPartitionsCreatesFourSurfaces)
        {
            gmsh::model::add("test");
            Geometry geo;

            std::vector<int> surfaces = geo.rectangleWithPartitions(4.0, 4.0, 0.0, 0.0, 0.1);

            // Should create 4 partitions (2x2 grid)
            EXPECT_EQ(surfaces.size(), 4u);

            // All surface tags should be positive
            for (int tag : surfaces)
            {
                EXPECT_GT(tag, 0);
            }
        }

        TEST_F(GeometryTest, RectangleWithPartitionsBoundingBox)
        {
            gmsh::model::add("test");
            Geometry geo;

            double length = 6.0, width = 4.0;
            double x = 1.0, y = 2.0;

            geo.rectangleWithPartitions(length, width, x, y, 0.1);
            BoundingBox bbox = geo.getBoundingBox();

            // OCC kernel has lower precision than GEO kernel, use 1e-6 tolerance
            EXPECT_NEAR(bbox[0], x, 1e-6);
            EXPECT_NEAR(bbox[1], y, 1e-6);
            EXPECT_NEAR(bbox[3], x + length, 1e-6);
            EXPECT_NEAR(bbox[4], y + width, 1e-6);
        }

        // =============================================================================
        // Multiple Geometry Tests
        // =============================================================================

        TEST_F(GeometryTest, MultipleRectangles)
        {
            gmsh::model::add("test");
            Geometry geo;

            int rect1 = geo.rectangle(1.0, 1.0, 0.0, 0.0, 0.1);
            int rect2 = geo.rectangle(1.0, 1.0, 2.0, 0.0, 0.1);

            EXPECT_GT(rect1, 0);
            EXPECT_GT(rect2, 0);
            EXPECT_NE(rect1, rect2);
        }

        TEST_F(GeometryTest, MixedGeometries)
        {
            gmsh::model::add("test");
            Geometry geo;

            int rect = geo.rectangle(2.0, 1.0, 0.0, 0.0, 0.1);
            int circ = geo.circle(0.5, 5.0, 0.5, 0.1);

            Point2D p1 = {7.0, 0.0};
            Point2D p2 = {8.0, 0.0};
            Point2D p3 = {7.5, 1.0};
            int tri = geo.triangle(p1, p2, p3, 0.1);

            EXPECT_GT(rect, 0);
            EXPECT_GT(circ, 0);
            EXPECT_GT(tri, 0);

            // All tags should be unique
            EXPECT_NE(rect, circ);
            EXPECT_NE(circ, tri);
            EXPECT_NE(rect, tri);
        }

        // =============================================================================
        // Mesh Size Tests
        // =============================================================================

        TEST_F(GeometryTest, DifferentMeshSizes)
        {
            gmsh::model::add("test");
            Geometry geo;

            // These should all succeed with different mesh sizes
            int s1 = geo.rectangle(1.0, 1.0, 0.0, 0.0, 0.01);
            EXPECT_GT(s1, 0);

            gmsh::clear();
            gmsh::model::add("test2");

            int s2 = geo.rectangle(1.0, 1.0, 0.0, 0.0, 1.0);
            EXPECT_GT(s2, 0);
        }

    } // namespace testing
} // namespace fvm
