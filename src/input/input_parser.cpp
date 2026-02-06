/**
 * @file input_parser.cpp
 * @brief Implementation of YAML parsing for mesh generation configuration.
 */

#include "input_parser.hpp"
#include "input_config.hpp"

#include <yaml-cpp/yaml.h>

#include <fstream>
#include <sstream>
#include <stdexcept>

namespace fvm
{

    namespace
    {

        // Helper to get a value with default
        template <typename T>
        T getWithDefault(const YAML::Node &node, const std::string &key, const T &defaultValue)
        {
            if (node[key])
            {
                return node[key].as<T>();
            }
            return defaultValue;
        }

        // Parse geometry configuration
        GeometryConfig parseGeometry(const YAML::Node &node)
        {
            GeometryConfig geo;

            geo.id = getWithDefault<std::string>(node, "id", "");
            geo.type = getWithDefault<std::string>(node, "type", "");

            if (node["params"])
            {
                const auto &params = node["params"];

                // Rectangle/ellipse common parameters
                geo.length = getWithDefault<double>(params, "length", 0.0);
                geo.width = getWithDefault<double>(params, "width", 0.0);
                geo.x = getWithDefault<double>(params, "x", 0.0);
                geo.y = getWithDefault<double>(params, "y", 0.0);

                // Circle
                geo.radius = getWithDefault<double>(params, "radius", 0.0);

                // Ellipse
                geo.r1 = getWithDefault<double>(params, "r1", 0.0);
                geo.r2 = getWithDefault<double>(params, "r2", 0.0);

                // Triangle vertices
                if (params["p1"] && params["p1"].IsSequence() && params["p1"].size() >= 2)
                {
                    geo.p1[0] = params["p1"][0].as<double>();
                    geo.p1[1] = params["p1"][1].as<double>();
                }
                if (params["p2"] && params["p2"].IsSequence() && params["p2"].size() >= 2)
                {
                    geo.p2[0] = params["p2"][0].as<double>();
                    geo.p2[1] = params["p2"][1].as<double>();
                }
                if (params["p3"] && params["p3"].IsSequence() && params["p3"].size() >= 2)
                {
                    geo.p3[0] = params["p3"][0].as<double>();
                    geo.p3[1] = params["p3"][1].as<double>();
                }

                // Polygon points
                if (params["points"] && params["points"].IsSequence())
                {
                    for (const auto &pt : params["points"])
                    {
                        if (pt.IsSequence() && pt.size() >= 2)
                        {
                            Point2D point = {pt[0].as<double>(), pt[1].as<double>()};
                            geo.points.push_back(point);
                        }
                    }
                }

                geo.convexHull = getWithDefault<bool>(params, "convexHull", false);
            }

            return geo;
        }

        // Parse boundary configuration
        BoundaryConfig parseBoundary(const YAML::Node &node)
        {
            BoundaryConfig bc;
            bc.name = getWithDefault<std::string>(node, "name", "");
            bc.expr = getWithDefault<std::string>(node, "expr", "");
            return bc;
        }

        /// Parse string to MeshType
        inline std::string stringToMeshType(const std::string &str)
        {
            if (str == "tri" || str == "triangle" || str == "triangles")
                return "tri";
            if (str == "quad" || str == "quads" || str == "quadrilateral")
                return "quad";
            if (str == "structured" || str == "structure")
                return "structured";
            return "tri"; // default
        }

        // Parse mesh configuration
        MeshConfig parseMesh(const YAML::Node &node)
        {
            MeshConfig mesh;

            std::string typeValue = getWithDefault<std::string>(node, "type", "triangles");

            mesh.meshType = stringToMeshType(typeValue);

            mesh.meshSize = getWithDefault<double>(node, "meshSize", 0.05);
            mesh.charLength = getWithDefault<double>(node, "charLength", 0.01);

            return mesh;
        }

        // Parse reorder configuration
        ReorderConfig parseReorder(const YAML::Node &node)
        {
            ReorderConfig reorder;

            reorder.cellStrategy = getWithDefault<std::string>(node, "cellStrategy", "");
            reorder.nodeStrategy = getWithDefault<std::string>(node, "nodeStrategy", "");

            return reorder;
        }

        // Parse partition configuration
        PartitionConfig parsePartition(const YAML::Node &node)
        {
            PartitionConfig partition;

            partition.enabled = getWithDefault<bool>(node, "enabled", false);
            partition.numParts = getWithDefault<int>(node, "numParts", 1);
            partition.method = getWithDefault<std::string>(node, "method", "metis");

            if (node["reorder"])
            {
                partition.reorder = parseReorder(node["reorder"]);
            }

            return partition;
        }

        // Parse output configuration
        OutputConfig parseOutput(const YAML::Node &node)
        {
            OutputConfig output;

            output.directory = getWithDefault<std::string>(node, "directory", "output");
            output.baseName = getWithDefault<std::string>(node, "baseName", "mesh");

            if (node["formats"] && node["formats"].IsSequence())
            {
                output.formats.clear();
                for (const auto &fmt : node["formats"])
                {
                    output.formats.push_back(fmt.as<std::string>());
                }
            }

            output.writeBoundaryInfo = getWithDefault<bool>(node, "writeBoundaryInfo", true);
            output.writePartitionMetadata = getWithDefault<bool>(node, "writePartitionMetadata", true);
            output.writeQualityReport = getWithDefault<bool>(node, "writeQualityReport", true);

            return output;
        }

        // Parse complete configuration from YAML node
        InputConfig parseConfig(const YAML::Node &root)
        {
            InputConfig config;

            // Project section
            if (root["project"])
            {
                const auto &project = root["project"];
                config.projectName = getWithDefault<std::string>(project, "name", "");
                config.projectDescription = getWithDefault<std::string>(project, "description", "");
            }

            // Geometries section
            if (root["geometries"] && root["geometries"].IsSequence())
            {
                for (const auto &geoNode : root["geometries"])
                {
                    config.geometries.push_back(parseGeometry(geoNode));
                }
            }
            else if (root["geometry"])
            {
                // Support single geometry for simple cases
                config.geometries.push_back(parseGeometry(root["geometry"]));
            }

            // Mesh section
            if (root["mesh"])
            {
                config.mesh = parseMesh(root["mesh"]);
            }

            // Boundaries section
            if (root["boundaries"] && root["boundaries"].IsSequence())
            {
                for (const auto &bcNode : root["boundaries"])
                {
                    config.boundaries.push_back(parseBoundary(bcNode));
                }
            }

            // Partition section
            if (root["partition"])
            {
                config.partition = parsePartition(root["partition"]);
            }

            // Reorder section
            if (root["reorder"])
            {
                config.reorder = parseReorder(root["reorder"]);
            }

            // Output section
            if (root["output"])
            {
                config.output = parseOutput(root["output"]);
            }

            return config;
        }

    } // anonymous namespace

    InputConfig parseYamlFile(const std::string &filepath)
    {
        try
        {
            YAML::Node root = YAML::LoadFile(filepath);
            return parseConfig(root);
        }
        catch (const YAML::Exception &e)
        {
            throw std::runtime_error("YAML parse error in '" + filepath + "': " + e.what());
        }
        catch (const std::exception &e)
        {
            throw std::runtime_error("Error reading '" + filepath + "': " + e.what());
        }
    }

    InputConfig parseYamlString(const std::string &yamlContent)
    {
        try
        {
            YAML::Node root = YAML::Load(yamlContent);
            return parseConfig(root);
        }
        catch (const YAML::Exception &e)
        {
            throw std::runtime_error(std::string("YAML parse error: ") + e.what());
        }
    }

} // namespace fvm
