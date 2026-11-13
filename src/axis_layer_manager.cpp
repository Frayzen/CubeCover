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
}

#include <algorithm>
#include <vector>

void AxisLayerManager::updateTraversalCuts(RegionId regionId) {
  const Axis nextAxis = (axis + 1) % DIMENSION;
  const Axis prevAxis = (axis + 2) % DIMENSION;

  Region &region = mesh.regions[regionId];

  const SliceId axisStart = region.startSlices(axis);
  const SliceId axisEnd = region.endSlices(axis);

  std::vector<std::array<bool, 2>> occupiedSlices(axisEnd - axisStart,
                                                  {false, false});

  auto processAxis = [&](Axis otherAxis, size_t occupidSliceId) {
    auto &manager = mesh.layerManagers[otherAxis];

    const RegionLayerId regionLayer = region.layerIds(otherAxis);

    for (RegionId adjacentId : manager.layerToRegions[regionLayer]) {
      if (adjacentId == regionId)
        continue;
      Region &adjacent = mesh.regions[adjacentId];
      bool isAdjacent = adjacent.layerIds[axis] == region.layerIds[axis];
      if (!isAdjacent)
        continue;

      const Axis remainingAxis = otherAxis == nextAxis ? prevAxis : nextAxis;
      SliceId face;
      if (adjacent.endSlices(remainingAxis) ==
          region.startSlices(remainingAxis)) {
        face = adjacent.endSlices(remainingAxis) - 1;
      } else if (adjacent.startSlices(remainingAxis) ==
                 region.endSlices(remainingAxis)) {
        face = adjacent.startSlices(remainingAxis);
      } else
        continue;

      const SliceId tStart = std::max(region.startSlices(otherAxis),
                                      adjacent.startSlices(otherAxis));

      const SliceId tEnd =
          std::min(region.endSlices(otherAxis), adjacent.endSlices(otherAxis));

      const SliceId sStart = std::max(axisStart, adjacent.startSlices(axis));
      const SliceId sEnd = std::min(axisEnd, adjacent.endSlices(axis));

      for (SliceId s = sStart; s < sEnd; ++s) {
        for (SliceId t = tStart; t < tEnd; ++t) {
          Eigen::Vector3i voxel = Eigen::Vector3i::Zero();

          voxel(axis) = s;
          voxel(remainingAxis) = face;
          voxel(otherAxis) = t;

          if (mesh.voxels(voxel.x(), voxel.y(), voxel.z()) != NO_VOXEL) {
            occupiedSlices[s - axisStart][occupidSliceId] = true;
            break;
          }
        }
      }
    }
  };

  processAxis(nextAxis, 0);
  processAxis(prevAxis, 1);

  TraversalCuts cuts;

  bool inRange = false;
  SliceId rangeStart = axisStart;

  for (SliceId s = axisStart; s < axisEnd; ++s) {
    if (occupiedSlices[s - axisStart][0] && occupiedSlices[s - axisStart][1]) {
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

bool AxisLayerManager::isCutAllowed(RegionId regionId, SliceId slice) const {
  Region &r = mesh.regions[regionId];
  if (layerToRegions[r.layerIds[axis]].size() < 3)
    return false;
  assert(r.containsSlice(axis, slice));
  for (const TraversalCut &cut : regionToChildrenTraversalCuts[regionId]) {
    if (cut.contains(slice))
      return true;
  }
  return false;
};
