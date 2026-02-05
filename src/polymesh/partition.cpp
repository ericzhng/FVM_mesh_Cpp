#include "partition.hpp"
#include "poly_mesh.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <stdexcept>

// Conditionally include METIS
#ifdef USE_METIS
#include <metis.h>
#endif

namespace fvm
{

    // =========================================================================
    // Utility Functions
    // =========================================================================

    std::vector<std::vector<Index>> getAdjacencyList(const PolyMesh &mesh)
    {
        std::vector<std::vector<Index>> adjacency(mesh.nCells);

        for (auto i = 0; i < mesh.nCells; ++i)
        {
            for (auto neighbor : mesh.cellNeighbors[i])
            {
                if (neighbor >= 0 && neighbor != i)
                {
                    adjacency[i].push_back(neighbor);
                }
            }
            // Remove duplicates and sort
            std::sort(adjacency[i].begin(), adjacency[i].end());
            adjacency[i].erase(
                std::unique(adjacency[i].begin(), adjacency[i].end()),
                adjacency[i].end());
        }

        return adjacency;
    }

    bool isMetisAvailable()
    {
#ifdef USE_METIS
        return true;
#else
        return false;
#endif
    }

    // =========================================================================
    // Main Partitioning Function
    // =========================================================================

    std::vector<Index> partitionMesh(
        const PolyMesh &mesh,
        Index nParts,
        const std::string &method,
        const std::vector<Real> &cellWeights)
    {
        if (nParts <= 1)
        {
            return std::vector<Index>(mesh.nCells, 0);
        }

        if (mesh.nCells == 0)
        {
            return {};
        }

        if (method == "metis")
        {
            return partitionWithMetis(mesh, nParts, cellWeights);
        }
        else if (method == "hierarchical")
        {
            return partitionWithHierarchical(mesh, nParts, cellWeights);
        }
        else
        {
            throw std::invalid_argument("Unknown partition method: " + method);
        }
    }

    // =========================================================================
    // METIS Partitioning
    // =========================================================================

    std::vector<Index> partitionWithMetis(
        const PolyMesh &mesh,
        Index nParts,
        const std::vector<Real> &cellWeights)
    {
#ifdef USE_METIS
        if (mesh.nCells == 0)
        {
            return {};
        }

        // Build adjacency list
        auto adjacency = getAdjacencyList(mesh);

        // Convert to METIS CSR format
        idx_t nvtxs = static_cast<idx_t>(mesh.nCells);
        std::vector<idx_t> xadj(nvtxs + 1);
        std::vector<idx_t> adjncy;

        xadj[0] = 0;
        for (idx_t i = 0; i < nvtxs; ++i)
        {
            xadj[i + 1] = xadj[i] + static_cast<idx_t>(adjacency[i].size());
            for (auto neighbor : adjacency[i])
            {
                adjncy.push_back(static_cast<idx_t>(neighbor));
            }
        }

        // Prepare vertex weights if provided
        std::vector<idx_t> vwgt;
        if (!cellWeights.empty())
        {
            vwgt.reserve(nvtxs);
            for (auto i = 0; i < mesh.nCells; ++i)
            {
                // METIS requires integer weights, scale appropriately
                vwgt.push_back(static_cast<idx_t>(cellWeights[i] * 1000 + 1));
            }
        }

        // METIS parameters
        idx_t ncon = 1;
        idx_t nparts = static_cast<idx_t>(nParts);
        idx_t objval;
        std::vector<idx_t> part(nvtxs);

        // Options
        idx_t options[METIS_NOPTIONS];
        METIS_SetDefaultOptions(options);
        options[METIS_OPTION_NUMBERING] = 0; // C-style numbering

        // Call METIS
        int status = METIS_PartGraphKway(
            &nvtxs,
            &ncon,
            xadj.data(),
            adjncy.data(),
            vwgt.empty() ? nullptr : vwgt.data(), // vertex weights
            nullptr,                              // vertex sizes
            nullptr,                              // edge weights
            &nparts,
            nullptr, // target partition weights
            nullptr, // ubvec
            options,
            &objval,
            part.data());

        if (status != METIS_OK)
        {
            throw std::runtime_error("METIS partitioning failed with status: " +
                                     std::to_string(status));
        }

        return std::vector<Index>(part.begin(), part.end());

#else
        // METIS not available, fall back to hierarchical
        std::cerr << "Warning: METIS not available, falling back to hierarchical partitioning\n";
        return partitionWithHierarchical(mesh, nParts, cellWeights);
#endif
    }

    // =========================================================================
    // Hierarchical Coordinate Bisection
    // =========================================================================

    std::vector<Index> partitionWithHierarchical(
        const PolyMesh &mesh,
        Index nParts,
        const std::vector<Real> &cellWeights)
    {
        if (mesh.nCells == 0)
        {
            return {};
        }

        // Check if nParts is power of two
        bool isPowerOfTwo = (nParts > 0) && ((nParts & (nParts - 1)) == 0);
        if (!isPowerOfTwo)
        {
            std::cerr << "Warning: Hierarchical method works best with power-of-two "
                      << "partitions. n_parts=" << nParts << " may result in uneven partitions.\n";
        }

        // Initialize weights (uniform if not provided)
        std::vector<Real> weights(mesh.nCells, 1.0);
        if (!cellWeights.empty() && cellWeights.size() == mesh.nCells)
        {
            weights = cellWeights;
        }

        // Initialize all cells to partition 0
        std::vector<Index> parts(mesh.nCells, 0);

        // Iteratively bisect partitions until we have nParts
        for (auto i = 1; i < nParts; ++i)
        {
            // Find the partition with the most cells to split
            std::vector<Index> partCounts(i + 1, 0);
            for (auto p : parts)
            {
                if (p < partCounts.size())
                {
                    partCounts[p]++;
                }
            }

            Index partToSplit = static_cast<Index>(
                std::max_element(partCounts.begin(), partCounts.end()) - partCounts.begin());

            // Get indices of cells in this partition
            std::vector<Index> idxsToSplit;
            for (auto j = 0; j < mesh.nCells; ++j)
            {
                if (parts[j] == partToSplit)
                {
                    idxsToSplit.push_back(j);
                }
            }

            if (idxsToSplit.empty())
                continue;

            // Determine axis to split (longest dimension of bounding box)
            Real minX = std::numeric_limits<Real>::max();
            Real maxX = std::numeric_limits<Real>::lowest();
            Real minY = std::numeric_limits<Real>::max();
            Real maxY = std::numeric_limits<Real>::lowest();
            Real minZ = std::numeric_limits<Real>::max();
            Real maxZ = std::numeric_limits<Real>::lowest();

            for (auto idx : idxsToSplit)
            {
                const auto &c = mesh.cellCentroids[idx];
                minX = std::min(minX, c[0]);
                maxX = std::max(maxX, c[0]);
                minY = std::min(minY, c[1]);
                maxY = std::max(maxY, c[1]);
                minZ = std::min(minZ, c[2]);
                maxZ = std::max(maxZ, c[2]);
            }

            Real rangeX = maxX - minX;
            Real rangeY = maxY - minY;
            Real rangeZ = maxZ - minZ;

            Index axis = 0; // Default to X
            if (rangeY > rangeX && rangeY >= rangeZ)
            {
                axis = 1;
            }
            else if (rangeZ > rangeX && rangeZ > rangeY)
            {
                axis = 2;
            }

            // Sort indices by coordinate along chosen axis
            std::sort(idxsToSplit.begin(), idxsToSplit.end(),
                      [&mesh, axis](std::size_t a, std::size_t b)
                      {
                          return mesh.cellCentroids[a][axis] < mesh.cellCentroids[b][axis];
                      });

            // Find weighted median split point
            std::vector<Real> sortedWeights;
            for (auto idx : idxsToSplit)
            {
                sortedWeights.push_back(weights[idx]);
            }

            Real totalWeight = std::accumulate(sortedWeights.begin(),
                                               sortedWeights.end(), 0.0);
            Real halfWeight = totalWeight / 2.0;

            std::size_t splitIdx = idxsToSplit.size() / 2;
            if (totalWeight > 0)
            {
                Real cumWeight = 0.0;
                for (auto k = 0; k < sortedWeights.size(); ++k)
                {
                    cumWeight += sortedWeights[k];
                    if (cumWeight >= halfWeight)
                    {
                        splitIdx = k;
                        break;
                    }
                }
            }

            // Handle edge cases
            if (splitIdx == 0)
                splitIdx = 1;
            if (splitIdx >= idxsToSplit.size())
                splitIdx = idxsToSplit.size() - 1;

            // Assign cells on the right side to new partition
            for (auto k = splitIdx; k < idxsToSplit.size(); ++k)
            {
                parts[idxsToSplit[k]] = i;
            }
        }

        return parts;
    }

    // =========================================================================
    // Partition Summary
    // =========================================================================

    void printPartitionSummary(const std::vector<Index> &parts)
    {
        if (parts.empty())
        {
            std::cout << "--- Partition Summary ---\n";
            std::cout << "No partitions found.\n";
            return;
        }

        Index nParts = *std::max_element(parts.begin(), parts.end()) + 1;

        // Count cells per partition
        std::vector<std::size_t> counts(nParts, 0);
        for (auto p : parts)
        {
            if (p >= 0 && p < nParts)
            {
                counts[p]++;
            }
        }

        std::cout << "--- Partition Summary ---\n";
        std::cout << "Number of partitions: " << nParts << "\n";

        std::size_t totalCells = parts.size();
        std::size_t minCells = *std::min_element(counts.begin(), counts.end());
        std::size_t maxCells = *std::max_element(counts.begin(), counts.end());
        Real avgCells = static_cast<Real>(totalCells) / nParts;

        for (auto p = 0; p < nParts; ++p)
        {
            Real pct = 100.0 * counts[p] / totalCells;
            std::cout << "  Partition " << std::setw(3) << p << ": "
                      << std::setw(8) << counts[p] << " cells ("
                      << std::fixed << std::setprecision(1) << pct << "%)\n";
        }

        std::cout << "\nBalance Statistics:\n";
        std::cout << "  Min cells per partition: " << minCells << "\n";
        std::cout << "  Max cells per partition: " << maxCells << "\n";
        std::cout << "  Avg cells per partition: " << std::fixed
                  << std::setprecision(1) << avgCells << "\n";

        Real imbalance = (maxCells > 0)
                             ? (static_cast<Real>(maxCells) / avgCells - 1.0) * 100.0
                             : 0.0;
        std::cout << "  Load imbalance: " << std::fixed << std::setprecision(1)
                  << imbalance << "%\n";
    }

} // namespace fvm
