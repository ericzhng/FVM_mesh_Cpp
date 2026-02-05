#define _USE_MATH_DEFINES
#include <cmath>

#include "mesh_quality.hpp"
#include "poly_mesh.hpp"

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <map>
#include <numeric>
#include <set>
#include <stdexcept>

namespace fvm
{

    // =========================================================================
    // Factory Method
    // =========================================================================

    MeshQuality MeshQuality::fromMesh(const PolyMesh &mesh)
    {
        if (!mesh.isAnalyzed())
        {
            throw std::runtime_error("Mesh must be analyzed before computing quality.");
        }

        MeshQuality quality;

        if (mesh.nCells == 0)
        {
            return quality;
        }

        quality.minMaxVolumeRatio = computeVolumeRatio(mesh);

        auto [skewness, aspectRatio] = computeGeometricMetrics(mesh);
        quality.cellSkewness = std::move(skewness);
        quality.cellAspectRatio = std::move(aspectRatio);

        quality.cellNonOrthogonality = computeNonOrthogonality(mesh);
        quality.connectivityIssues = checkConnectivity(mesh);

        return quality;
    }

    // =========================================================================
    // Output Methods
    // =========================================================================

    void MeshQuality::printSummary() const
    {
        std::cout << "\n"
                  << std::setw(40) << "--- Quality Metrics ---" << "\n\n";

        std::cout << std::fixed << std::setprecision(4);
        std::cout << "  " << std::left << std::setw(30)
                  << "Min/Max Volume Ratio:" << minMaxVolumeRatio << "\n";

        // Helper lambda for statistics
        auto printStats = [](const std::string &name, const std::vector<Real> &data)
        {
            if (data.empty())
                return;

            std::vector<Real> validData;
            for (Real v : data)
            {
                if (std::isfinite(v))
                    validData.push_back(v);
            }
            if (validData.empty())
                return;

            Real minVal = *std::min_element(validData.begin(), validData.end());
            Real maxVal = *std::max_element(validData.begin(), validData.end());
            Real avgVal = std::accumulate(validData.begin(), validData.end(), 0.0) /
                          static_cast<Real>(validData.size());

            std::cout << "  " << std::left << std::setw(25) << name
                      << std::right << std::setw(12) << minVal
                      << std::setw(12) << maxVal
                      << std::setw(12) << avgVal << "\n";
        };

        std::cout << "\n  " << std::left << std::setw(25) << "Metric"
                  << std::right << std::setw(12) << "Min"
                  << std::setw(12) << "Max"
                  << std::setw(12) << "Average" << "\n";
        std::cout << "  " << std::string(24, '-') << " "
                  << std::string(12, '-') << " "
                  << std::string(12, '-') << " "
                  << std::string(12, '-') << "\n";

        printStats("Skewness", cellSkewness);
        printStats("Aspect Ratio", cellAspectRatio);
        printStats("Non-Orthogonality (deg)", cellNonOrthogonality);

        if (!connectivityIssues.empty())
        {
            std::cout << "\n  Connectivity Issues:\n";
            for (const auto &issue : connectivityIssues)
            {
                std::cout << "    - " << issue << "\n";
            }
        }
        else
        {
            std::cout << "\n  Connectivity: OK (no issues found)\n";
        }
    }

    // =========================================================================
    // Internal Computation Methods
    // =========================================================================

    Real MeshQuality::computeVolumeRatio(const PolyMesh &mesh)
    {
        if (mesh.cellVolumes.empty())
            return 0.0;

        Real minVol = *std::min_element(mesh.cellVolumes.begin(), mesh.cellVolumes.end());
        Real maxVol = *std::max_element(mesh.cellVolumes.begin(), mesh.cellVolumes.end());

        return (maxVol > GEOMETRY_TOLERANCE) ? (minVol / maxVol) : 0.0;
    }

    std::pair<std::vector<Real>, std::vector<Real>>
    MeshQuality::computeGeometricMetrics(const PolyMesh &mesh)
    {
        Index nCells = mesh.nCells;
        std::vector<Real> skewness(nCells, 0.0);
        std::vector<Real> aspectRatio(nCells, 1.0);

        // Currently only supports 2D meshes with triangles and quads
        if (mesh.dimension != 2)
        {
            return {skewness, aspectRatio};
        }

        for (auto i = 0; i < nCells; ++i)
        {
            const auto &conn = mesh.cellNodeConnectivity[i];
            auto numNodes = conn.size();

            // Extract 2D node coordinates
            std::vector<std::array<Real, 2>> nodes2D;
            nodes2D.reserve(numNodes);
            for (auto nodeIdx : conn)
            {
                nodes2D.push_back({mesh.nodeCoords[nodeIdx][0],
                                   mesh.nodeCoords[nodeIdx][1]});
            }

            if (numNodes == 3)
            {
                auto [s, ar] = computeTriangleMetrics(nodes2D);
                skewness[i] = s;
                aspectRatio[i] = ar;
            }
            else if (numNodes == 4)
            {
                auto [s, ar] = computeQuadMetrics(nodes2D);
                skewness[i] = s;
                aspectRatio[i] = ar;
            }
        }

        return {skewness, aspectRatio};
    }

    std::pair<Real, Real>
    MeshQuality::computeTriangleMetrics(const std::vector<std::array<Real, 2>> &nodes)
    {
        // Compute edge lengths
        std::vector<Real> lengths(3);
        for (int i = 0; i < 3; ++i)
        {
            int j = (i + 1) % 3;
            Real dx = nodes[j][0] - nodes[i][0];
            Real dy = nodes[j][1] - nodes[i][1];
            lengths[i] = std::sqrt(dx * dx + dy * dy);
        }

        Real minLen = *std::min_element(lengths.begin(), lengths.end());
        if (minLen < GEOMETRY_TOLERANCE)
        {
            return {1.0, std::numeric_limits<Real>::infinity()};
        }

        // Compute angles at each vertex using law of cosines
        std::vector<Real> angles(3);
        for (int i = 0; i < 3; ++i)
        {
            int prev = (i + 2) % 3; // Previous vertex
            int next = (i + 1) % 3; // Next vertex

            // Vectors from vertex i to neighbors
            Real v1x = nodes[prev][0] - nodes[i][0];
            Real v1y = nodes[prev][1] - nodes[i][1];
            Real v2x = nodes[next][0] - nodes[i][0];
            Real v2y = nodes[next][1] - nodes[i][1];

            Real len1 = lengths[prev]; // Edge from i to prev
            Real len2 = lengths[i];    // Edge from i to next

            Real dot = v1x * v2x + v1y * v2y;
            Real cosAngle = std::clamp(dot / (len1 * len2), -1.0, 1.0);
            angles[i] = std::acos(cosAngle) * 180.0 / M_PI;
        }

        // Skewness: max deviation from ideal angle (60 degrees)
        Real maxDeviation = 0.0;
        for (Real angle : angles)
        {
            maxDeviation = std::max(maxDeviation, std::abs(angle - IDEAL_TRIANGLE_ANGLE));
        }
        Real skewness = maxDeviation / IDEAL_TRIANGLE_ANGLE;

        // Aspect ratio: longest edge / shortest edge
        Real maxLen = *std::max_element(lengths.begin(), lengths.end());
        Real aspectRatio = maxLen / minLen;

        return {skewness, aspectRatio};
    }

    std::pair<Real, Real>
    MeshQuality::computeQuadMetrics(const std::vector<std::array<Real, 2>> &nodes)
    {
        // Compute edge lengths
        std::vector<Real> lengths(4);
        for (int i = 0; i < 4; ++i)
        {
            int j = (i + 1) % 4;
            Real dx = nodes[j][0] - nodes[i][0];
            Real dy = nodes[j][1] - nodes[i][1];
            lengths[i] = std::sqrt(dx * dx + dy * dy);
        }

        Real minLen = *std::min_element(lengths.begin(), lengths.end());
        if (minLen < GEOMETRY_TOLERANCE)
        {
            return {1.0, std::numeric_limits<Real>::infinity()};
        }

        // Compute angles at each vertex
        std::vector<Real> angles(4);
        for (int i = 0; i < 4; ++i)
        {
            int prev = (i + 3) % 4;
            int next = (i + 1) % 4;

            Real v1x = nodes[prev][0] - nodes[i][0];
            Real v1y = nodes[prev][1] - nodes[i][1];
            Real v2x = nodes[next][0] - nodes[i][0];
            Real v2y = nodes[next][1] - nodes[i][1];

            Real len1 = std::sqrt(v1x * v1x + v1y * v1y);
            Real len2 = std::sqrt(v2x * v2x + v2y * v2y);

            if (len1 < GEOMETRY_TOLERANCE || len2 < GEOMETRY_TOLERANCE)
            {
                return {1.0, std::numeric_limits<Real>::infinity()};
            }

            Real dot = v1x * v2x + v1y * v2y;
            Real cosAngle = std::clamp(dot / (len1 * len2), -1.0, 1.0);
            angles[i] = std::acos(cosAngle) * 180.0 / M_PI;
        }

        // Skewness: max deviation from ideal angle (90 degrees)
        Real maxDeviation = 0.0;
        for (Real angle : angles)
        {
            maxDeviation = std::max(maxDeviation, std::abs(angle - IDEAL_QUAD_ANGLE));
        }
        Real skewness = maxDeviation / IDEAL_QUAD_ANGLE;

        // Aspect ratio: longest edge / shortest edge
        Real maxLen = *std::max_element(lengths.begin(), lengths.end());
        Real aspectRatio = maxLen / minLen;

        return {skewness, aspectRatio};
    }

    std::vector<Real> MeshQuality::computeNonOrthogonality(const PolyMesh &mesh)
    {
        std::vector<Real> nonOrtho(mesh.nCells, 0.0);

        for (auto ci = 0; ci < mesh.nCells; ++ci)
        {
            Real maxNonOrtho = 0.0;

            for (auto fi = 0; fi < mesh.cellFaceNodes[ci].size(); ++fi)
            {
                if (fi >= mesh.cellFaceMidpoints[ci].size())
                    continue;

                // Vector from cell centroid to face midpoint
                Real vecX = mesh.cellFaceMidpoints[ci][fi][0] - mesh.cellCentroids[ci][0];
                Real vecY = mesh.cellFaceMidpoints[ci][fi][1] - mesh.cellCentroids[ci][1];
                Real vecZ = mesh.cellFaceMidpoints[ci][fi][2] - mesh.cellCentroids[ci][2];

                Real normVec = std::sqrt(vecX * vecX + vecY * vecY + vecZ * vecZ);

                // Face normal
                Real normX = mesh.cellFaceNormals[ci][fi][0];
                Real normY = mesh.cellFaceNormals[ci][fi][1];
                Real normZ = mesh.cellFaceNormals[ci][fi][2];

                Real normNormal = std::sqrt(normX * normX + normY * normY + normZ * normZ);

                if (normVec > GEOMETRY_TOLERANCE && normNormal > GEOMETRY_TOLERANCE)
                {
                    Real dot = vecX * normX + vecY * normY + vecZ * normZ;
                    Real cosAngle = std::clamp(dot / (normVec * normNormal), -1.0, 1.0);
                    Real angleDeg = std::acos(cosAngle) * 180.0 / M_PI;
                    maxNonOrtho = std::max(maxNonOrtho, angleDeg);
                }
            }

            nonOrtho[ci] = maxNonOrtho;
        }

        return nonOrtho;
    }

    std::vector<std::string> MeshQuality::checkConnectivity(const PolyMesh &mesh)
    {
        std::vector<std::string> issues;

        // Check for unreferenced nodes
        std::set<Index> referencedNodes;
        for (const auto &conn : mesh.cellNodeConnectivity)
        {
            for (auto nodeIdx : conn)
            {
                referencedNodes.insert(nodeIdx);
            }
        }

        if (referencedNodes.size() < mesh.nNodes)
        {
            auto unreferenced = mesh.nNodes - referencedNodes.size();
            issues.push_back("Found " + std::to_string(unreferenced) + " unreferenced nodes.");
        }

        // Check for duplicate cells
        std::set<std::vector<Index>> uniqueCells;
        for (const auto &conn : mesh.cellNodeConnectivity)
        {
            std::vector<Index> sortedConn = conn;
            std::sort(sortedConn.begin(), sortedConn.end());
            uniqueCells.insert(sortedConn);
        }

        if (uniqueCells.size() < mesh.nCells)
        {
            issues.push_back("Found duplicate cells.");
        }

        // Check for non-manifold faces (faces shared by more than 2 cells)
        std::map<std::vector<Index>, std::vector<Index>> faceMap;
        for (auto ci = 0; ci < mesh.nCells; ++ci)
        {
            for (const auto &face : mesh.cellFaceNodes[ci])
            {
                std::vector<Index> key = face;
                std::sort(key.begin(), key.end());
                faceMap[key].push_back(ci);
            }
        }

        for (const auto &[key, cells] : faceMap)
        {
            if (cells.size() > 2)
            {
                std::string faceStr = "(";
                for (auto i = 0; i < key.size(); ++i)
                {
                    if (i > 0)
                        faceStr += ", ";
                    faceStr += std::to_string(key[i]);
                }
                faceStr += ")";
                issues.push_back("Non-manifold face " + faceStr +
                                 " shared by " + std::to_string(cells.size()) + " cells.");
            }
        }

        return issues;
    }

} // namespace fvm
