#include "mesh_generator.hpp"
#include <gmsh.h>
#include <algorithm>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <set>
#include <stdexcept>
#include <unordered_set>
#include "vtkio/cell_types.hpp"

namespace fvm
{

    MeshGenerator::MeshGenerator(const std::vector<int> &surfaceTags,
                                 const std::string &outputDir)
        : surfaceTags_(surfaceTags), outputDir_(outputDir)
    {
        ensureOutputDir();
    }

    MeshGenerator::MeshGenerator(int surfaceTag, const std::string &outputDir)
        : surfaceTags_({surfaceTag}), outputDir_(outputDir)
    {
        ensureOutputDir();
    }

    void MeshGenerator::ensureOutputDir()
    {
        std::filesystem::create_directories(outputDir_);
    }

    void MeshGenerator::generate(const std::map<int, MeshParams> &meshParams,
                                 const std::string &filename)
    {
        // Apply mesh parameters for each surface
        for (auto surfaceTag : surfaceTags_)
        {
            auto it = meshParams.find(surfaceTag);
            if (it != meshParams.end())
            {
                applyMeshParameters(surfaceTag, it->second);
            }
        }

        // Set up physical groups
        setupPhysicalGroups();

        // Generate 2D mesh
        gmsh::model::mesh::generate(2);

        // Save mesh to Gmsh format
        saveMesh(filename);

        // Extract mesh data for further processing
        extractMeshData();

        std::cout << "Mesh generation complete." << std::endl;
        std::cout << "  Nodes: " << meshData_.nodes.size() << std::endl;
        std::cout << "  Elements: " << meshData_.elements.size() << std::endl;
    }

    void MeshGenerator::applyMeshParameters(int surfaceTag, const MeshParams &params)
    {
        if (params.meshType == "structured")
        {
            setStructuredMesh(surfaceTag, params.charLength);
        }
        else if (params.meshType == "quad")
        {
            gmsh::model::mesh::setRecombine(2, surfaceTag);
        }
        else
        {
            // do nothing for triangle
        }
    }

    void MeshGenerator::setStructuredMesh(int surfaceTag, double charLength)
    {
        gmsh::model::mesh::setTransfiniteSurface(surfaceTag);
        gmsh::model::mesh::setRecombine(2, surfaceTag);

        // Get boundary curves
        std::vector<std::pair<int, int>> boundaryCurves;
        gmsh::model::getBoundary({{2, surfaceTag}}, boundaryCurves, false, false, false);

        if (boundaryCurves.size() != 4)
        {
            throw std::runtime_error(
                "Structured mesh is only supported for geometries with 4 boundary curves.");
        }

        // Get bounding box for calculating number of divisions
        double minX, minY, minZ, maxX, maxY, maxZ;
        gmsh::model::getBoundingBox(2, surfaceTag, minX, minY, minZ, maxX, maxY, maxZ);

        double dx = maxX - minX;
        double dy = maxY - minY;
        auto nx = static_cast<int>(dx / charLength);
        auto ny = static_cast<int>(dy / charLength);

        // Classify curves and set transfinite
        auto [hCurves, vCurves] = classifyBoundaryCurves(boundaryCurves);

        for (int curveTag : hCurves)
        {
            gmsh::model::mesh::setTransfiniteCurve(curveTag, nx + 1);
        }
        for (int curveTag : vCurves)
        {
            gmsh::model::mesh::setTransfiniteCurve(curveTag, ny + 1);
        }
    }

    std::pair<std::vector<int>, std::vector<int>>
    MeshGenerator::classifyBoundaryCurves(
        const std::vector<std::pair<int, int>> &boundaryCurves) const
    {
        std::vector<int> hCurves, vCurves;

        for (const auto &dimTag : boundaryCurves)
        {
            int curveTag = dimTag.second;

            // Get curve endpoints
            std::vector<std::pair<int, int>> pointTags;
            gmsh::model::getBoundary({{1, curveTag}}, pointTags, false, false, false);

            if (pointTags.size() < 2)
                continue;

            int pStartTag = pointTags[0].second;
            int pEndTag = pointTags[1].second;

            std::vector<double> coordStart, coordEnd;
            std::vector<double> paramCoord; // unused but required
            gmsh::model::getValue(0, pStartTag, paramCoord, coordStart);
            gmsh::model::getValue(0, pEndTag, paramCoord, coordEnd);

            // Classify by orientation
            if (std::abs(coordStart[1] - coordEnd[1]) < 1e-6)
            {
                hCurves.push_back(curveTag);
            }
            else
            {
                vCurves.push_back(curveTag);
            }
        }

        return {hCurves, vCurves};
    }

    void MeshGenerator::setupPhysicalGroups()
    {
        // Collect all boundary curves
        std::set<int> allBoundaryCurves;
        for (auto surfaceTag : surfaceTags_)
        {
            std::vector<std::pair<int, int>> boundary;
            gmsh::model::getBoundary({{2, surfaceTag}}, boundary, false, false, false);
            for (const auto &dimTag : boundary)
            {
                allBoundaryCurves.insert(dimTag.second);
            }
        }

        // Find curves already in physical groups
        std::set<int> curvesInGroups;
        std::vector<std::pair<int, int>> physicalGroups;
        gmsh::model::getPhysicalGroups(physicalGroups, 1);

        for (const auto &[dim, tag] : physicalGroups)
        {
            std::vector<int> entities;
            gmsh::model::getEntitiesForPhysicalGroup(dim, tag, entities);
            for (auto e : entities)
            {
                curvesInGroups.insert(e);
            }
        }

        // Add untagged curves to "unnamed" group
        std::vector<int> untaggedCurves;
        for (auto curve : allBoundaryCurves)
        {
            if (curvesInGroups.find(curve) == curvesInGroups.end())
            {
                untaggedCurves.push_back(curve);
            }
        }
        if (!untaggedCurves.empty())
        {
            gmsh::model::addPhysicalGroup(1, untaggedCurves, -1, "unnamed");
        }

        // Handle 2D physical groups (surfaces)
        std::set<int> surfacesInGroups;
        gmsh::model::getPhysicalGroups(physicalGroups, 2);

        for (const auto &[dim, tag] : physicalGroups)
        {
            std::vector<int> entities;
            gmsh::model::getEntitiesForPhysicalGroup(dim, tag, entities);
            for (auto e : entities)
            {
                surfacesInGroups.insert(e);
            }
        }

        // Add untagged surfaces to "fluid" group
        std::vector<int> untaggedSurfaces;
        for (auto surface : surfaceTags_)
        {
            if (surfacesInGroups.find(surface) == surfacesInGroups.end())
            {
                untaggedSurfaces.push_back(surface);
            }
        }
        if (!untaggedSurfaces.empty())
        {
            gmsh::model::addPhysicalGroup(2, untaggedSurfaces, -1, "fluid");
        }
    }

    void MeshGenerator::saveMesh(const std::string &filename)
    {
        std::string mshFile = outputDir_ + "/" + filename;
        gmsh::write(mshFile);
        std::cout << "Mesh saved to: " << mshFile << std::endl;
    }

    void MeshGenerator::extractMeshData()
    {
        extractNodes();
        extractCells();
        extractPhysicalGroups();
        extractFaces();
    }

    void MeshGenerator::extractNodes()
    {
        std::vector<std::size_t> nodeTags;
        std::vector<double> nodeCoords;
        std::vector<double> parametricCoords;

        gmsh::model::mesh::getNodes(nodeTags, nodeCoords, parametricCoords);

        meshData_.nodes.clear();
        nodeIds_.clear();
        meshData_.nodes.reserve(nodeTags.size());
        nodeIds_.reserve(nodeTags.size());

        for (auto i = 0; i < nodeTags.size(); ++i)
        {
            nodeIds_.push_back(nodeTags[i]);
            meshData_.nodes.push_back({nodeCoords[3 * i],
                                       nodeCoords[3 * i + 1],
                                       nodeCoords[3 * i + 2]});
        }
    }

    void MeshGenerator::extractCells()
    {
        // Create node tag to index map
        std::unordered_map<std::size_t, std::size_t> nodeMap;
        for (auto i = 0; i < nodeIds_.size(); ++i)
        {
            nodeMap[nodeIds_[i]] = i;
        }

        meshData_.elements.clear();
        meshData_.elementTypes.clear();

        for (auto surfaceTag : surfaceTags_)
        {
            std::vector<int> elemTypes;
            std::vector<std::vector<std::size_t>> elemTags;
            std::vector<std::vector<std::size_t>> elemNodeTags;

            gmsh::model::mesh::getElements(elemTypes, elemTags, elemNodeTags, 2, surfaceTag);

            for (auto i = 0; i < elemTypes.size(); ++i)
            {
                auto elemType = elemTypes[i];

                // Get element properties
                std::string name;
                int dim, order, numNodes, numPrimaryNodes;
                std::vector<double> localNodeCoords;
                gmsh::model::mesh::getElementProperties(
                    elemType, name, dim, order, numNodes, localNodeCoords, numPrimaryNodes);

                auto numElements = elemTags[i].size();
                const auto &allNodeTags = elemNodeTags[i];

                for (auto j = 0; j < numElements; ++j)
                {
                    CellConnectivity cell;
                    cell.reserve(numNodes);

                    for (auto k = 0; k < numNodes; ++k)
                    {
                        auto nodeTag = allNodeTags[j * numNodes + k];
                        cell.push_back(nodeMap[nodeTag]);
                    }

                    meshData_.elements.push_back(std::move(cell));
                    meshData_.elementTypes.push_back(getVTKCellType(numNodes));
                }
            }
        }
    }

    void MeshGenerator::extractPhysicalGroups()
    {
        meshData_.nodeSets.clear();
        meshData_.elementSets.clear();
        meshData_.faceSets.clear();

        // Create node tag to index map
        std::unordered_map<std::size_t, std::size_t> nodeMap;
        for (auto i = 0; i < nodeIds_.size(); ++i)
        {
            nodeMap[nodeIds_[i]] = i;
        }

        std::vector<std::pair<int, int>> physicalGroups;
        gmsh::model::getPhysicalGroups(physicalGroups);

        for (const auto &[dim, tag] : physicalGroups)
        {
            std::string name;
            gmsh::model::getPhysicalName(dim, tag, name);
            if (name.empty())
            {
                name = "group_" + std::to_string(tag);
            }

            std::vector<int> entities;
            gmsh::model::getEntitiesForPhysicalGroup(dim, tag, entities);

            if (dim == 1)
            {
                // Boundary curves -> face sets (edges in 2D)
                std::vector<FaceNodes> faces;
                for (auto entityTag : entities)
                {
                    std::vector<int> elemTypes;
                    std::vector<std::vector<std::size_t>> elemTags;
                    std::vector<std::vector<std::size_t>> elemNodeTags;

                    gmsh::model::mesh::getElements(elemTypes, elemTags, elemNodeTags, 1, entityTag);

                    for (auto i = 0; i < elemTypes.size(); ++i)
                    {
                        const auto &nodeTags = elemNodeTags[i];
                        // Line elements have 2 nodes
                        for (auto j = 0; j + 1 < nodeTags.size(); j += 2)
                        {
                            FaceNodes face = {static_cast<Index>(nodeMap[nodeTags[j]]), static_cast<Index>(nodeMap[nodeTags[j + 1]])};
                            faces.push_back(face);
                        }
                    }
                }
                meshData_.faceSets[name] = std::move(faces);
            }
            else if (dim == 2)
            {
                // Surface groups -> element sets
                std::vector<Index> elementIndices;
                // Note: For simplicity, we store entity tags; proper implementation
                // would map surface elements to their indices
                for (auto entityTag : entities)
                {
                    elementIndices.push_back(static_cast<Index>(entityTag));
                }
                meshData_.elementSets[name] = std::move(elementIndices);
            }
        }
    }

    void MeshGenerator::extractFaces()
    {
        // Boundary faces are already extracted in extractPhysicalGroups() as faceSets.
        // This function can be extended to compute internal faces if needed for FVM.
    }

} // namespace fvm
