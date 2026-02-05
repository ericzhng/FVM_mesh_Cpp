#pragma once

/**
 * @file fvm_types.hpp
 * @brief Common type definitions shared across FVM modules.
 *
 * This header contains fundamental types used by both meshgen and polymesh
 * modules, allowing them to remain independent libraries.
 *
 * Numeric precision is controlled by macros in fvm_export.hpp:
 * - FVM_USE_SINGLE_PRECISION: Use float instead of double
 * - FVM_USE_32BIT_INT: Use 32-bit integers instead of 64-bit
 */

#include "fvm_export.hpp"
#include "vtkio/cell_types.hpp" // VTK cell type definitions
#include <array>
#include <string>
#include <unordered_map>
#include <vector>

namespace fvm
{
    // =============================================================================
    // Basic Geometric Types (using configurable Real type)
    // =============================================================================

    /// 2D point/vector type
    using Point2D = std::array<Real, 2>;

    /// 3D point/vector type
    using Point3D = std::array<Real, 3>;

    /// Bounding box: {min_x, min_y, min_z, max_x, max_y, max_z}
    using BoundingBox = std::array<Real, 6>;

    /// Cell connectivity (variable number of nodes)
    using CellConnectivity = std::vector<Id>;

    /// A face defined by its node indices (variable size: 3 for tri, 4 for quad, etc.)
    using FaceNodes = std::vector<Id>;

    // =============================================================================
    // Mesh Interchange Data Structure
    // =============================================================================

    /// Minimal mesh data structure for format interchange
    struct MeshInfo
    {
        // === Geometry ===
        std::vector<Point3D> nodes;

        // === Topology ===
        std::vector<CellConnectivity> elements;
        std::vector<Sid> elementTypes; // VTK cell types (-1 if unknown)

        // === Named Sets ===
        std::unordered_map<std::string, std::vector<Id>> nodeSets;
        std::unordered_map<std::string, std::vector<Id>> elementSets;
        std::unordered_map<std::string, std::vector<FaceNodes>> faceSets; // boundary groups go here
    };

    // =============================================================================
    // Element and Physical Group Definitions
    // =============================================================================

    /**
     * @brief Properties of an element type (e.g., triangle, quad, tetrahedron).
     */
    struct ElementTypeProperties
    {
        std::string name;
        int numNodes;
    };

    /**
     * @brief Information about a physical group (boundary or volume region).
     */
    struct PhysicalGroup
    {
        Sid dimension;
        Sid tag;
        std::string name;
        std::vector<Sid> entities;
    };

    // =============================================================================
    // Mesh Type Enumeration
    // =============================================================================

    /// Mesh type enumeration
    enum class MeshType
    {
        Triangles, // Unstructured triangular mesh
        Quads,     // Unstructured quadrilateral mesh
        Structured // Structured quadrilateral mesh
    };

    /// Convert MeshType to string
    inline std::string meshTypeToString(MeshType type)
    {
        switch (type)
        {
        case MeshType::Triangles:
            return "tri";
        case MeshType::Quads:
            return "quads";
        case MeshType::Structured:
            return "structured";
        default:
            return "unknown";
        }
    }

    /// Parse string to MeshType
    inline MeshType stringToMeshType(const std::string &str)
    {
        if (str == "tri" || str == "triangles")
            return MeshType::Triangles;
        if (str == "quads" || str == "quad")
            return MeshType::Quads;
        if (str == "structured")
            return MeshType::Structured;
        return MeshType::Triangles; // default
    }

    /// Mesh parameters for a surface
    struct MeshParams
    {
        MeshType meshType = MeshType::Triangles;
        Real charLength = 0.1;
    };

} // namespace fvm
