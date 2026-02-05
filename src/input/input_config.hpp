#pragma once

/**
 * @file input_config.hpp
 * @brief Configuration data structures for unified mesh generation input.
 *
 * This header defines the structures used to represent a complete mesh
 * generation configuration, including geometry, mesh parameters, boundaries,
 * partitioning, and output settings.
 */

#include "common/fvm_export.hpp"
#include "common/fvm_types.hpp"

#include <string>
#include <vector>

namespace fvm
{

    // =============================================================================
    // Geometry Configuration
    // =============================================================================

    /**
     * @brief Configuration for a single geometry shape.
     *
     * This union-like struct holds parameters for all supported geometry types.
     * Only the fields relevant to the specified type need to be set.
     */
    struct FVM_API GeometryConfig
    {
        std::string id;   ///< Optional identifier for the geometry
        std::string type; ///< Geometry type: rectangle, circle, triangle, ellipse, polygon

        // Rectangle parameters
        double length = 0.0; ///< Length (x-direction) for rectangle
        double width = 0.0;  ///< Width (y-direction) for rectangle

        // Circle parameters
        double radius = 0.0; ///< Radius for circle

        // Ellipse parameters
        double r1 = 0.0; ///< X-axis radius for ellipse
        double r2 = 0.0; ///< Y-axis radius for ellipse

        // Common position parameters
        double x = 0.0; ///< X position (center for circle/ellipse, corner for rectangle)
        double y = 0.0; ///< Y position

        // Triangle parameters
        Point2D p1 = {0.0, 0.0}; ///< First vertex for triangle
        Point2D p2 = {0.0, 0.0}; ///< Second vertex for triangle
        Point2D p3 = {0.0, 0.0}; ///< Third vertex for triangle

        // Polygon parameters
        std::vector<Point2D> points; ///< Vertices for polygon
        bool convexHull = false;     ///< Whether to compute convex hull for polygon
    };

    // =============================================================================
    // Boundary Configuration
    // =============================================================================

    /**
     * @brief Configuration for a boundary condition.
     *
     * Boundaries are identified using mathematical expressions evaluated at
     * edge midpoints. The expression should evaluate to true (non-zero) for
     * edges that belong to this boundary.
     *
     * Supported operators: ==, <, >, <=, >=, +, -, *, /, ^, ||, &&
     * Supported variables: X, Y (edge midpoint coordinates)
     * Supported functions: sqrt(), abs()
     *
     * Examples:
     *   - "X == 0" (left boundary)
     *   - "Y == 1.0" (top boundary)
     *   - "X^2 + Y^2 <= 0.25" (inside circle of radius 0.5)
     */
    struct FVM_API BoundaryConfig
    {
        std::string name; ///< Name of the boundary (e.g., "inlet", "outlet", "wall")
        std::string expr; ///< Mathematical expression for edge identification
    };

    // =============================================================================
    // Mesh Configuration
    // =============================================================================

    /**
     * @brief Configuration for mesh generation parameters.
     */
    struct FVM_API MeshConfig
    {
        std::string meshType = ""; ///< Type of mesh elements
        double meshSize = 0.05;    ///< Initial mesh size for geometry creation
        double charLength = 0.01;  ///< Characteristic element length
    };

    // =============================================================================
    // Partition Configuration
    // =============================================================================

    /**
     * @brief Configuration for mesh partitioning.
     */
    struct FVM_API PartitionConfig
    {
        bool enabled = false;         ///< Whether to partition the mesh
        int numParts = 1;             ///< Number of partitions
        std::string method = "metis"; ///< Partitioning method: "metis" or "hierarchical"
    };

    // =============================================================================
    // Reorder Configuration
    // =============================================================================

    /**
     * @brief Configuration for cell and node reordering within partitions.
     *
     * Supported cell strategies: rcm, gps, sloan, spectral, spatial_x, spatial_y, random
     * Supported node strategies: rcm, sequential, reverse, spatial_x, spatial_y, random
     */
    struct FVM_API ReorderConfig
    {
        std::string cellStrategy; ///< Cell reordering strategy (empty for no reordering)
        std::string nodeStrategy; ///< Node reordering strategy (empty for no reordering)
    };

    // =============================================================================
    // Output Configuration
    // =============================================================================

    /**
     * @brief Configuration for output files and formats.
     */
    struct FVM_API OutputConfig
    {
        std::string directory = "output";           ///< Output directory path
        std::string baseName = "mesh";              ///< Base name for output files
        std::vector<std::string> formats = {"vtu"}; ///< Export formats: vtu, vtk, msh, openfoam
        bool writeBoundaryInfo = true;              ///< Write boundary information file
        bool writePartitionMetadata = true;         ///< Write partition metadata (JSON)
    };

    // =============================================================================
    // Complete Input Configuration
    // =============================================================================

    /**
     * @brief Complete configuration for mesh generation.
     *
     * This structure holds all settings needed to generate a mesh from
     * geometry definitions through partitioning and export.
     */
    struct FVM_API InputConfig
    {
        // Project information
        std::string projectName;        ///< Project name
        std::string projectDescription; ///< Optional description

        // Geometry definitions (multiple geometries form a union)
        std::vector<GeometryConfig> geometries;

        // Mesh settings
        MeshConfig mesh;

        // Boundary definitions
        std::vector<BoundaryConfig> boundaries;

        // Partitioning settings
        PartitionConfig partition;

        // Reordering settings
        ReorderConfig reorder;

        // Output settings
        OutputConfig output;

        /**
         * @brief Validate the configuration.
         * @param errorMessage Output parameter for error description
         * @return true if configuration is valid, false otherwise
         */
        bool validate(std::string &errorMessage) const;
    };

} // namespace fvm
