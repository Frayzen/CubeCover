#include "axis_layer_manager.hpp"
#include "ids.hpp"
#include "voxel_mesh.hpp"

RegionLayerId AxisLayerManager::prepareSplit(RegionLayerId layerId) {
  RegionLayerId newLayer = childGraph.createNode();
  assert(parentGraph.createNode() == newLayer);
  auto layerChildren = childGraph[layerId];
  auto layerParents = parentGraph[layerId];
  for (RegionId parent : layerParents) {
    auto &parentChildren = childGraph[parent];
    for (size_t i = 0; i < parentChildren.size(); i++)
      if (parentChildren[i] == layerId)
        parentChildren[i] = newLayer;
  }
  childGraph[newLayer] = {layerId};
  parentGraph[layerId] = {newLayer};
  parentGraph[newLayer] = layerParents;

  layerToRegions.push_back({});
  return newLayer;
}

void AxisLayerManager::splitRegionFromAnotherAxis(RegionId splittedRegionId,
                                                  RegionId regionId) {
  RegionLayerId layerId = mesh.regions[splittedRegionId].layerIds[axis];

  layerToRegions[layerId].push_back(regionId);
  regionToChildrenTraversalCuts.push_back({});
  // TODO smarter way splitting cuts depending on the axis
  updateTraversalCuts(regionId);
  updateTraversalCuts(splittedRegionId);
}

#include <algorithm>
#include <vector>

void AxisLayerManager::updateTraversalCuts(RegionId regionId) {
  const Axis nextAxis = (axis + 1) % DIMENSION;
  const Axis prevAxis = (axis + 2) % DIMENSION;

  Region &region = mesh.regions[regionId];

  const SliceId axisStart = region.startSlices(axis);
  const SliceId axisEnd = region.endSlices(axis);

  std::vector<bool> occupiedSlices(axisEnd - axisStart, false);

  auto processAxis = [&](Axis otherAxis) {
    auto &manager = mesh.layerManagers[otherAxis];

    const RegionLayerId parentLayer = region.layerIds(otherAxis);

    // Children in the other axis layer manager.
    for (RegionLayerId adjacentLayer : manager.childGraph[parentLayer]) {

      for (RegionId adjacentId : manager.layerToRegions[adjacentLayer]) {
        Region &adjacent = mesh.regions[adjacentId];

        // Determine whether the adjacent region shares a face
        // perpendicular to otherAxis.
        SliceId face;

        if (adjacent.endSlices(otherAxis) == region.startSlices(otherAxis)) {
          face = adjacent.endSlices(otherAxis) - 1;
        } else if (adjacent.startSlices(otherAxis) ==
                   region.endSlices(otherAxis)) {
          face = adjacent.startSlices(otherAxis);
        } else {
          continue;
        }

        const Axis remainingAxis = otherAxis == nextAxis ? prevAxis : nextAxis;

        const SliceId tStart = std::max(region.startSlices(remainingAxis),
                                        adjacent.startSlices(remainingAxis));

        const SliceId tEnd = std::min(region.endSlices(remainingAxis),
                                      adjacent.endSlices(remainingAxis));

        const SliceId sStart = std::max(axisStart, adjacent.startSlices(axis));

        const SliceId sEnd = std::min(axisEnd, adjacent.endSlices(axis));

        for (SliceId s = sStart; s < sEnd; ++s) {
          for (SliceId t = tStart; t < tEnd; ++t) {
            Eigen::Vector3i voxel = Eigen::Vector3i::Zero();

            voxel(axis) = s;
            voxel(otherAxis) = face;
            voxel(remainingAxis) = t;

            if (mesh.voxels(voxel.x(), voxel.y(), voxel.z()) != NO_VOXEL) {
              occupiedSlices[s - axisStart] = true;
              break;
            }
          }
        }
      }
    }
  };

  // Process both other axis layer managers.
  processAxis(nextAxis);
  processAxis(prevAxis);

  // Merge consecutive occupied slices.
  TraversalCuts cuts;

  bool inRange = false;
  SliceId rangeStart = axisStart;

  for (SliceId s = axisStart; s < axisEnd; ++s) {
    if (occupiedSlices[s - axisStart]) {
      if (!inRange) {
        rangeStart = s;
        inRange = true;
      }
    } else if (inRange) {
      cuts.push_back({rangeStart, s});
      inRange = false;
    }
  }

  if (inRange)
    cuts.push_back({rangeStart, axisEnd});

  regionToChildrenTraversalCuts[regionId] = std::move(cuts);
}

void AxisLayerManager::concretizeSplit(RegionLayerId layerId,
                                       std::vector<RegionId> &regions) {
  layerToRegions[layerId] = regions;
  for (size_t i = 0; i < regions.size(); i++)
    regionToChildrenTraversalCuts.push_back({});

  assert(childGraph[layerId].size() == 1);
  RegionLayerId childLayer = childGraph[layerId][0];
  for (RegionId rId : layerToRegions[layerId])
    updateTraversalCuts(rId);
  for (RegionId rId : layerToRegions[childLayer])
    updateTraversalCuts(rId);
}

void AxisLayerManager::initialize(Region &initRegion) {
  parentGraph.createNode();
  childGraph.createNode();
  TraversalCut cc = {
      .from = initRegion.startSlices[axis],
      .to = initRegion.endSlices[axis],
  };
  regionToChildrenTraversalCuts.push_back({cc});
  layerToRegions.push_back({0});
}

bool AxisLayerManager::isCutAllowed(RegionLayerId id, SliceId slice) const {
  if (layerToRegions[id].size() < 4)
    return false;
  for (RegionId regionId : layerToRegions[id]) {
    for (auto cut : regionToChildrenTraversalCuts[regionId])
      if (!cut.contains(slice))
        return false;
  }
  return true;
};
