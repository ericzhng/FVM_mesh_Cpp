/**
 * @file input_config.cpp
 * @brief Implementation of configuration validation.
 */

#include "input_config.hpp"

#include <set>
#include <sstream>

namespace fvm
{

    bool InputConfig::validate(std::string &errorMessage) const
    {
        std::ostringstream errors;
        bool valid = true;

        // Validate project name
        if (projectName.empty())
        {
            errors << "- project.name is required\n";
            valid = false;
        }

        // Validate geometries
        if (geometries.empty())
        {
            errors << "- At least one geometry is required\n";
            valid = false;
        }

        static const std::set<std::string> validGeometryTypes = {
            "rectangle", "circle", "triangle", "ellipse", "polygon"};

        for (std::size_t i = 0; i < geometries.size(); ++i)
        {
            const auto &geo = geometries[i];

            // Check geometry type
            if (validGeometryTypes.find(geo.type) == validGeometryTypes.end())
            {
                errors << "- geometries[" << i << "]: invalid type '" << geo.type << "'\n";
                valid = false;
                continue;
            }

            // Type-specific validation
            if (geo.type == "rectangle")
            {
                if (geo.length <= 0)
                {
                    errors << "- geometries[" << i << "]: rectangle requires positive length\n";
                    valid = false;
                }
                if (geo.width <= 0)
                {
                    errors << "- geometries[" << i << "]: rectangle requires positive width\n";
                    valid = false;
                }
            }
            else if (geo.type == "circle")
            {
                if (geo.radius <= 0)
                {
                    errors << "- geometries[" << i << "]: circle requires positive radius\n";
                    valid = false;
                }
            }
            else if (geo.type == "ellipse")
            {
                if (geo.r1 <= 0 || geo.r2 <= 0)
                {
                    errors << "- geometries[" << i << "]: ellipse requires positive r1 and r2\n";
                    valid = false;
                }
            }
            else if (geo.type == "polygon")
            {
                if (geo.points.size() < 3)
                {
                    errors << "- geometries[" << i << "]: polygon requires at least 3 points\n";
                    valid = false;
                }
            }
            // Triangle validation: check that points are not collinear
            // (simplified check - just ensure they're not all the same)
            else if (geo.type == "triangle")
            {
                if (geo.p1 == geo.p2 || geo.p2 == geo.p3 || geo.p1 == geo.p3)
                {
                    errors << "- geometries[" << i << "]: triangle vertices must be distinct\n";
                    valid = false;
                }
            }
        }

        // Validate mesh parameters
        if (mesh.meshSize <= 0)
        {
            errors << "- mesh.meshSize must be positive\n";
            valid = false;
        }
        if (mesh.charLength <= 0)
        {
            errors << "- mesh.charLength must be positive\n";
            valid = false;
        }

        // Validate boundaries
        for (std::size_t i = 0; i < boundaries.size(); ++i)
        {
            const auto &bc = boundaries[i];
            if (bc.name.empty())
            {
                errors << "- boundaries[" << i << "]: name is required\n";
                valid = false;
            }
            if (bc.expr.empty())
            {
                errors << "- boundaries[" << i << "]: expr is required\n";
                valid = false;
            }
        }

        // Validate partition settings
        if (partition.enabled)
        {
            if (partition.numParts < 1)
            {
                errors << "- partition.numParts must be at least 1\n";
                valid = false;
            }
            if (partition.method != "metis" && partition.method != "hierarchical")
            {
                errors << "- partition.method must be 'metis' or 'hierarchical'\n";
                valid = false;
            }
        }

        // Validate reorder strategies
        static const std::set<std::string> validCellStrategies = {
            "", "rcm", "gps", "sloan", "spectral", "spatial_x", "spatial_y", "random"};
        static const std::set<std::string> validNodeStrategies = {
            "", "rcm", "sequential", "reverse", "spatial_x", "spatial_y", "random"};

        if (validCellStrategies.find(reorder.cellStrategy) == validCellStrategies.end())
        {
            errors << "- reorder.cellStrategy '" << reorder.cellStrategy << "' is not valid\n";
            valid = false;
        }
        if (validNodeStrategies.find(reorder.nodeStrategy) == validNodeStrategies.end())
        {
            errors << "- reorder.nodeStrategy '" << reorder.nodeStrategy << "' is not valid\n";
            valid = false;
        }

        // Validate output settings
        if (output.directory.empty())
        {
            errors << "- output.directory is required\n";
            valid = false;
        }
        if (output.baseName.empty())
        {
            errors << "- output.baseName is required\n";
            valid = false;
        }

        static const std::set<std::string> validFormats = {
            "vtu", "vtk", "msh", "openfoam"};
        for (const auto &fmt : output.formats)
        {
            if (validFormats.find(fmt) == validFormats.end())
            {
                errors << "- output.formats: unknown format '" << fmt << "'\n";
                valid = false;
            }
        }

        errorMessage = errors.str();
        return valid;
    }

} // namespace fvm
