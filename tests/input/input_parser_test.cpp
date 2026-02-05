#include <gtest/gtest.h>
#include "input/input_parser.hpp"
#include "input/input_config.hpp"

namespace fvm
{
    namespace testing
    {

        // =============================================================================
        // YAML Parsing Tests
        // =============================================================================

        TEST(InputParserYaml, ParseMinimalConfig)
        {
            std::string yaml = R"(
project:
  name: "Test Project"
geometries:
  - type: rectangle
    params:
      length: 1.0
      width: 1.0
)";
            InputConfig config = parseYamlString(yaml);
            EXPECT_EQ(config.projectName, "Test Project");
            EXPECT_EQ(config.geometries.size(), 1u);
            EXPECT_EQ(config.geometries[0].type, "rectangle");
        }

        TEST(InputParserYaml, ParseRectangleGeometry)
        {
            std::string yaml = R"(
project:
  name: "Test"
geometries:
  - id: "main_rect"
    type: rectangle
    params:
      length: 10.0
      width: 5.0
      x: 1.0
      y: 2.0
)";
            InputConfig config = parseYamlString(yaml);

            ASSERT_EQ(config.geometries.size(), 1u);
            const auto &geo = config.geometries[0];
            EXPECT_EQ(geo.id, "main_rect");
            EXPECT_EQ(geo.type, "rectangle");
            EXPECT_DOUBLE_EQ(geo.length, 10.0);
            EXPECT_DOUBLE_EQ(geo.width, 5.0);
            EXPECT_DOUBLE_EQ(geo.x, 1.0);
            EXPECT_DOUBLE_EQ(geo.y, 2.0);
        }

        TEST(InputParserYaml, ParseCircleGeometry)
        {
            std::string yaml = R"(
project:
  name: "Test"
geometries:
  - type: circle
    params:
      radius: 5.0
      x: 10.0
      y: 20.0
)";
            InputConfig config = parseYamlString(yaml);

            ASSERT_EQ(config.geometries.size(), 1u);
            const auto &geo = config.geometries[0];
            EXPECT_EQ(geo.type, "circle");
            EXPECT_DOUBLE_EQ(geo.radius, 5.0);
            EXPECT_DOUBLE_EQ(geo.x, 10.0);
            EXPECT_DOUBLE_EQ(geo.y, 20.0);
        }

        TEST(InputParserYaml, ParseEllipseGeometry)
        {
            std::string yaml = R"(
project:
  name: "Test"
geometries:
  - type: ellipse
    params:
      r1: 3.0
      r2: 2.0
      x: 1.0
      y: 1.0
)";
            InputConfig config = parseYamlString(yaml);

            ASSERT_EQ(config.geometries.size(), 1u);
            const auto &geo = config.geometries[0];
            EXPECT_EQ(geo.type, "ellipse");
            EXPECT_DOUBLE_EQ(geo.r1, 3.0);
            EXPECT_DOUBLE_EQ(geo.r2, 2.0);
        }

        TEST(InputParserYaml, ParseTriangleGeometry)
        {
            std::string yaml = R"(
project:
  name: "Test"
geometries:
  - type: triangle
    params:
      p1: [0.0, 0.0]
      p2: [1.0, 0.0]
      p3: [0.5, 1.0]
)";
            InputConfig config = parseYamlString(yaml);

            ASSERT_EQ(config.geometries.size(), 1u);
            const auto &geo = config.geometries[0];
            EXPECT_EQ(geo.type, "triangle");
            EXPECT_DOUBLE_EQ(geo.p1[0], 0.0);
            EXPECT_DOUBLE_EQ(geo.p1[1], 0.0);
            EXPECT_DOUBLE_EQ(geo.p2[0], 1.0);
            EXPECT_DOUBLE_EQ(geo.p2[1], 0.0);
            EXPECT_DOUBLE_EQ(geo.p3[0], 0.5);
            EXPECT_DOUBLE_EQ(geo.p3[1], 1.0);
        }

        TEST(InputParserYaml, ParsePolygonGeometry)
        {
            std::string yaml = R"(
project:
  name: "Test"
geometries:
  - type: polygon
    params:
      points:
        - [0.0, 0.0]
        - [1.0, 0.0]
        - [1.0, 1.0]
        - [0.0, 1.0]
      convexHull: true
)";
            InputConfig config = parseYamlString(yaml);

            ASSERT_EQ(config.geometries.size(), 1u);
            const auto &geo = config.geometries[0];
            EXPECT_EQ(geo.type, "polygon");
            EXPECT_EQ(geo.points.size(), 4u);
            EXPECT_TRUE(geo.convexHull);
        }

        TEST(InputParserYaml, ParseMeshSettings)
        {
            std::string yaml = R"(
project:
  name: "Test"
geometries:
  - type: rectangle
    params:
      length: 1.0
      width: 1.0
mesh:
  type: quads
  meshSize: 0.1
  charLength: 0.05
)";
            InputConfig config = parseYamlString(yaml);

            EXPECT_EQ(config.mesh.meshType, MeshType::Quads);
            EXPECT_DOUBLE_EQ(config.mesh.meshSize, 0.1);
            EXPECT_DOUBLE_EQ(config.mesh.charLength, 0.05);
        }

        TEST(InputParserYaml, ParseTriangleMeshType)
        {
            std::string yaml = R"(
project:
  name: "Test"
geometries:
  - type: rectangle
    params:
      length: 1.0
      width: 1.0
mesh:
  type: triangles
)";
            InputConfig config = parseYamlString(yaml);
            EXPECT_EQ(config.mesh.meshType, MeshType::Triangles);
        }

        TEST(InputParserYaml, ParseBoundaries)
        {
            std::string yaml = R"(
project:
  name: "Test"
geometries:
  - type: rectangle
    params:
      length: 1.0
      width: 1.0
boundaries:
  - name: inlet
    expr: "X == 0"
  - name: outlet
    expr: "X == 1"
  - name: wall
    expr: "Y == 0 || Y == 1"
)";
            InputConfig config = parseYamlString(yaml);

            ASSERT_EQ(config.boundaries.size(), 3u);
            EXPECT_EQ(config.boundaries[0].name, "inlet");
            EXPECT_EQ(config.boundaries[0].expr, "X == 0");
            EXPECT_EQ(config.boundaries[1].name, "outlet");
            EXPECT_EQ(config.boundaries[1].expr, "X == 1");
            EXPECT_EQ(config.boundaries[2].name, "wall");
            EXPECT_EQ(config.boundaries[2].expr, "Y == 0 || Y == 1");
        }

        TEST(InputParserYaml, ParsePartitionSettings)
        {
            std::string yaml = R"(
project:
  name: "Test"
geometries:
  - type: rectangle
    params:
      length: 1.0
      width: 1.0
partition:
  enabled: true
  numParts: 4
  method: metis
)";
            InputConfig config = parseYamlString(yaml);

            EXPECT_TRUE(config.partition.enabled);
            EXPECT_EQ(config.partition.numParts, 4);
            EXPECT_EQ(config.partition.method, "metis");
        }

        TEST(InputParserYaml, ParseReorderSettings)
        {
            std::string yaml = R"(
project:
  name: "Test"
geometries:
  - type: rectangle
    params:
      length: 1.0
      width: 1.0
reorder:
  cellStrategy: rcm
  nodeStrategy: sequential
)";
            InputConfig config = parseYamlString(yaml);

            EXPECT_EQ(config.reorder.cellStrategy, "rcm");
            EXPECT_EQ(config.reorder.nodeStrategy, "sequential");
        }

        TEST(InputParserYaml, ParseOutputSettings)
        {
            std::string yaml = R"(
project:
  name: "Test"
geometries:
  - type: rectangle
    params:
      length: 1.0
      width: 1.0
output:
  directory: "results"
  baseName: "my_mesh"
  formats:
    - vtu
    - vtk
  writeBoundaryInfo: false
  writePartitionMetadata: true
)";
            InputConfig config = parseYamlString(yaml);

            EXPECT_EQ(config.output.directory, "results");
            EXPECT_EQ(config.output.baseName, "my_mesh");
            ASSERT_EQ(config.output.formats.size(), 2u);
            EXPECT_EQ(config.output.formats[0], "vtu");
            EXPECT_EQ(config.output.formats[1], "vtk");
            EXPECT_FALSE(config.output.writeBoundaryInfo);
            EXPECT_TRUE(config.output.writePartitionMetadata);
        }

        TEST(InputParserYaml, ParseMultipleGeometries)
        {
            std::string yaml = R"(
project:
  name: "Test"
geometries:
  - type: rectangle
    params:
      length: 2.0
      width: 1.0
  - type: circle
    params:
      radius: 0.5
      x: 1.0
      y: 0.5
)";
            InputConfig config = parseYamlString(yaml);

            ASSERT_EQ(config.geometries.size(), 2u);
            EXPECT_EQ(config.geometries[0].type, "rectangle");
            EXPECT_EQ(config.geometries[1].type, "circle");
        }

        TEST(InputParserYaml, ParseSingleGeometryShorthand)
        {
            std::string yaml = R"(
project:
  name: "Test"
geometry:
  type: rectangle
  params:
    length: 1.0
    width: 1.0
)";
            InputConfig config = parseYamlString(yaml);
            ASSERT_EQ(config.geometries.size(), 1u);
            EXPECT_EQ(config.geometries[0].type, "rectangle");
        }

        TEST(InputParserYaml, DefaultValues)
        {
            std::string yaml = R"(
project:
  name: "Test"
geometries:
  - type: rectangle
    params:
      length: 1.0
      width: 1.0
)";
            InputConfig config = parseYamlString(yaml);

            // Check default mesh values
            EXPECT_EQ(config.mesh.meshType, MeshType::Triangles);
            EXPECT_DOUBLE_EQ(config.mesh.meshSize, 0.05);
            EXPECT_DOUBLE_EQ(config.mesh.charLength, 0.01);

            // Check default partition values
            EXPECT_FALSE(config.partition.enabled);
            EXPECT_EQ(config.partition.numParts, 1);
            EXPECT_EQ(config.partition.method, "metis");

            // Check default output values
            EXPECT_EQ(config.output.directory, "output");
            EXPECT_EQ(config.output.baseName, "mesh");
            EXPECT_TRUE(config.output.writeBoundaryInfo);
            EXPECT_TRUE(config.output.writePartitionMetadata);
        }

        TEST(InputParserYaml, InvalidYamlThrows)
        {
            std::string yaml = "this: is: invalid: yaml:";
            EXPECT_THROW(parseYamlString(yaml), std::runtime_error);
        }

        // =============================================================================
        // Configuration Validation Tests
        // =============================================================================

        TEST(InputConfigValidation, ValidMinimalConfig)
        {
            InputConfig config;
            config.projectName = "Test";
            config.geometries.push_back({});
            config.geometries[0].type = "rectangle";
            config.geometries[0].length = 1.0;
            config.geometries[0].width = 1.0;

            std::string error;
            EXPECT_TRUE(config.validate(error));
        }

        TEST(InputConfigValidation, MissingProjectName)
        {
            InputConfig config;
            config.geometries.push_back({});
            config.geometries[0].type = "rectangle";
            config.geometries[0].length = 1.0;
            config.geometries[0].width = 1.0;

            std::string error;
            EXPECT_FALSE(config.validate(error));
            EXPECT_NE(error.find("project.name"), std::string::npos);
        }

        TEST(InputConfigValidation, NoGeometries)
        {
            InputConfig config;
            config.projectName = "Test";

            std::string error;
            EXPECT_FALSE(config.validate(error));
            EXPECT_NE(error.find("geometry"), std::string::npos);
        }

        TEST(InputConfigValidation, InvalidGeometryType)
        {
            InputConfig config;
            config.projectName = "Test";
            config.geometries.push_back({});
            config.geometries[0].type = "invalid_type";

            std::string error;
            EXPECT_FALSE(config.validate(error));
            EXPECT_NE(error.find("invalid type"), std::string::npos);
        }

        TEST(InputConfigValidation, RectangleZeroLength)
        {
            InputConfig config;
            config.projectName = "Test";
            config.geometries.push_back({});
            config.geometries[0].type = "rectangle";
            config.geometries[0].length = 0.0;
            config.geometries[0].width = 1.0;

            std::string error;
            EXPECT_FALSE(config.validate(error));
            EXPECT_NE(error.find("positive length"), std::string::npos);
        }

        TEST(InputConfigValidation, RectangleZeroWidth)
        {
            InputConfig config;
            config.projectName = "Test";
            config.geometries.push_back({});
            config.geometries[0].type = "rectangle";
            config.geometries[0].length = 1.0;
            config.geometries[0].width = 0.0;

            std::string error;
            EXPECT_FALSE(config.validate(error));
            EXPECT_NE(error.find("positive width"), std::string::npos);
        }

        TEST(InputConfigValidation, CircleZeroRadius)
        {
            InputConfig config;
            config.projectName = "Test";
            config.geometries.push_back({});
            config.geometries[0].type = "circle";
            config.geometries[0].radius = 0.0;

            std::string error;
            EXPECT_FALSE(config.validate(error));
            EXPECT_NE(error.find("positive radius"), std::string::npos);
        }

        TEST(InputConfigValidation, EllipseZeroR1)
        {
            InputConfig config;
            config.projectName = "Test";
            config.geometries.push_back({});
            config.geometries[0].type = "ellipse";
            config.geometries[0].r1 = 0.0;
            config.geometries[0].r2 = 1.0;

            std::string error;
            EXPECT_FALSE(config.validate(error));
            EXPECT_NE(error.find("positive r1 and r2"), std::string::npos);
        }

        TEST(InputConfigValidation, PolygonTooFewPoints)
        {
            InputConfig config;
            config.projectName = "Test";
            config.geometries.push_back({});
            config.geometries[0].type = "polygon";
            config.geometries[0].points = {{0, 0}, {1, 0}}; // Only 2 points

            std::string error;
            EXPECT_FALSE(config.validate(error));
            EXPECT_NE(error.find("at least 3 points"), std::string::npos);
        }

        TEST(InputConfigValidation, TriangleDuplicateVertices)
        {
            InputConfig config;
            config.projectName = "Test";
            config.geometries.push_back({});
            config.geometries[0].type = "triangle";
            config.geometries[0].p1 = {0.0, 0.0};
            config.geometries[0].p2 = {0.0, 0.0}; // Same as p1
            config.geometries[0].p3 = {1.0, 1.0};

            std::string error;
            EXPECT_FALSE(config.validate(error));
            EXPECT_NE(error.find("distinct"), std::string::npos);
        }

        TEST(InputConfigValidation, ZeroMeshSize)
        {
            InputConfig config;
            config.projectName = "Test";
            config.geometries.push_back({});
            config.geometries[0].type = "rectangle";
            config.geometries[0].length = 1.0;
            config.geometries[0].width = 1.0;
            config.mesh.meshSize = 0.0;

            std::string error;
            EXPECT_FALSE(config.validate(error));
            EXPECT_NE(error.find("meshSize"), std::string::npos);
        }

        TEST(InputConfigValidation, ZeroCharLength)
        {
            InputConfig config;
            config.projectName = "Test";
            config.geometries.push_back({});
            config.geometries[0].type = "rectangle";
            config.geometries[0].length = 1.0;
            config.geometries[0].width = 1.0;
            config.mesh.charLength = 0.0;

            std::string error;
            EXPECT_FALSE(config.validate(error));
            EXPECT_NE(error.find("charLength"), std::string::npos);
        }

        TEST(InputConfigValidation, BoundaryMissingName)
        {
            InputConfig config;
            config.projectName = "Test";
            config.geometries.push_back({});
            config.geometries[0].type = "rectangle";
            config.geometries[0].length = 1.0;
            config.geometries[0].width = 1.0;
            config.boundaries.push_back({});
            config.boundaries[0].expr = "X == 0";

            std::string error;
            EXPECT_FALSE(config.validate(error));
            EXPECT_NE(error.find("name is required"), std::string::npos);
        }

        TEST(InputConfigValidation, BoundaryMissingExpr)
        {
            InputConfig config;
            config.projectName = "Test";
            config.geometries.push_back({});
            config.geometries[0].type = "rectangle";
            config.geometries[0].length = 1.0;
            config.geometries[0].width = 1.0;
            config.boundaries.push_back({});
            config.boundaries[0].name = "inlet";

            std::string error;
            EXPECT_FALSE(config.validate(error));
            EXPECT_NE(error.find("expr is required"), std::string::npos);
        }

        TEST(InputConfigValidation, PartitionZeroParts)
        {
            InputConfig config;
            config.projectName = "Test";
            config.geometries.push_back({});
            config.geometries[0].type = "rectangle";
            config.geometries[0].length = 1.0;
            config.geometries[0].width = 1.0;
            config.partition.enabled = true;
            config.partition.numParts = 0;

            std::string error;
            EXPECT_FALSE(config.validate(error));
            EXPECT_NE(error.find("at least 1"), std::string::npos);
        }

        TEST(InputConfigValidation, PartitionInvalidMethod)
        {
            InputConfig config;
            config.projectName = "Test";
            config.geometries.push_back({});
            config.geometries[0].type = "rectangle";
            config.geometries[0].length = 1.0;
            config.geometries[0].width = 1.0;
            config.partition.enabled = true;
            config.partition.numParts = 4;
            config.partition.method = "invalid_method";

            std::string error;
            EXPECT_FALSE(config.validate(error));
            EXPECT_NE(error.find("method"), std::string::npos);
        }

        TEST(InputConfigValidation, InvalidCellStrategy)
        {
            InputConfig config;
            config.projectName = "Test";
            config.geometries.push_back({});
            config.geometries[0].type = "rectangle";
            config.geometries[0].length = 1.0;
            config.geometries[0].width = 1.0;
            config.reorder.cellStrategy = "invalid_strategy";

            std::string error;
            EXPECT_FALSE(config.validate(error));
            EXPECT_NE(error.find("cellStrategy"), std::string::npos);
        }

        TEST(InputConfigValidation, InvalidNodeStrategy)
        {
            InputConfig config;
            config.projectName = "Test";
            config.geometries.push_back({});
            config.geometries[0].type = "rectangle";
            config.geometries[0].length = 1.0;
            config.geometries[0].width = 1.0;
            config.reorder.nodeStrategy = "invalid_strategy";

            std::string error;
            EXPECT_FALSE(config.validate(error));
            EXPECT_NE(error.find("nodeStrategy"), std::string::npos);
        }

        TEST(InputConfigValidation, EmptyOutputDirectory)
        {
            InputConfig config;
            config.projectName = "Test";
            config.geometries.push_back({});
            config.geometries[0].type = "rectangle";
            config.geometries[0].length = 1.0;
            config.geometries[0].width = 1.0;
            config.output.directory = "";

            std::string error;
            EXPECT_FALSE(config.validate(error));
            EXPECT_NE(error.find("directory"), std::string::npos);
        }

        TEST(InputConfigValidation, EmptyBaseName)
        {
            InputConfig config;
            config.projectName = "Test";
            config.geometries.push_back({});
            config.geometries[0].type = "rectangle";
            config.geometries[0].length = 1.0;
            config.geometries[0].width = 1.0;
            config.output.baseName = "";

            std::string error;
            EXPECT_FALSE(config.validate(error));
            EXPECT_NE(error.find("baseName"), std::string::npos);
        }

        TEST(InputConfigValidation, InvalidOutputFormat)
        {
            InputConfig config;
            config.projectName = "Test";
            config.geometries.push_back({});
            config.geometries[0].type = "rectangle";
            config.geometries[0].length = 1.0;
            config.geometries[0].width = 1.0;
            config.output.formats = {"vtu", "invalid_format"};

            std::string error;
            EXPECT_FALSE(config.validate(error));
            EXPECT_NE(error.find("unknown format"), std::string::npos);
        }

        TEST(InputConfigValidation, ValidCellStrategies)
        {
            std::vector<std::string> validStrategies = {"", "rcm", "gps", "sloan", "spectral", "spatial_x", "spatial_y", "random"};

            for (const auto &strategy : validStrategies)
            {
                InputConfig config;
                config.projectName = "Test";
                config.geometries.push_back({});
                config.geometries[0].type = "rectangle";
                config.geometries[0].length = 1.0;
                config.geometries[0].width = 1.0;
                config.reorder.cellStrategy = strategy;

                std::string error;
                EXPECT_TRUE(config.validate(error)) << "Strategy '" << strategy << "' should be valid";
            }
        }

        TEST(InputConfigValidation, ValidNodeStrategies)
        {
            std::vector<std::string> validStrategies = {"", "rcm", "sequential", "reverse", "spatial_x", "spatial_y", "random"};

            for (const auto &strategy : validStrategies)
            {
                InputConfig config;
                config.projectName = "Test";
                config.geometries.push_back({});
                config.geometries[0].type = "rectangle";
                config.geometries[0].length = 1.0;
                config.geometries[0].width = 1.0;
                config.reorder.nodeStrategy = strategy;

                std::string error;
                EXPECT_TRUE(config.validate(error)) << "Strategy '" << strategy << "' should be valid";
            }
        }

        TEST(InputConfigValidation, ValidOutputFormats)
        {
            std::vector<std::string> validFormats = {"vtu", "vtk", "msh", "openfoam"};

            for (const auto &format : validFormats)
            {
                InputConfig config;
                config.projectName = "Test";
                config.geometries.push_back({});
                config.geometries[0].type = "rectangle";
                config.geometries[0].length = 1.0;
                config.geometries[0].width = 1.0;
                config.output.formats = {format};

                std::string error;
                EXPECT_TRUE(config.validate(error)) << "Format '" << format << "' should be valid";
            }
        }

    } // namespace testing
} // namespace fvm
