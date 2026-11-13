#include "axis_layer_manager.hpp"
#include "ids.hpp"
#include "voxel_mesh.hpp"
#include <algorithm>
#include <vector>

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

void AxisLayerManager::updateTraversalCuts(RegionId regionId) {
  const Axis nextAxis = (axis + 1) % DIMENSION;
  const Axis prevAxis = (axis + 2) % DIMENSION;
  const Axis axisList[2] = {nextAxis, prevAxis};

  Region &region = mesh.regions[regionId];

  const SliceId axisStart = region.startSlices(axis);
  const SliceId axisEnd = region.endSlices(axis);

  std::vector<std::array<bool, 2>> occupiedSlices(axisEnd - axisStart,
                                                  {false, false});

  auto processAxis = [&](size_t axisListId) {
    const Axis otherAxis = axisList[axisListId];
    const Axis remainingAxis = otherAxis == nextAxis ? prevAxis : nextAxis;

    const SliceId oAxisEnd = region.endSlices(otherAxis);
    const SliceId oAxisStart = region.startSlices(otherAxis);

    for (RegionId adjacentId : mesh.getAdjacentRegions(regionId, otherAxis)) {
      Region &adjacent = mesh.regions[adjacentId];
      const SliceId vStart = std::max(region.startSlices(remainingAxis),
                                      adjacent.startSlices(remainingAxis));
      const SliceId vEnd = std::min(region.endSlices(remainingAxis),
                                    adjacent.endSlices(remainingAxis));
      const SliceId uStart = std::max(axisStart, adjacent.startSlices(axis));
      const SliceId uEnd = std::min(axisEnd, adjacent.endSlices(axis));
      const SliceId face = adjacent.startSlices(otherAxis) == oAxisEnd
                               ? oAxisEnd
                               : oAxisStart - 1;
      for (SliceId u = uStart; u < uEnd; ++u) {
        Eigen::Vector3i voxel = Eigen::Vector3i::Zero();
        voxel(axis) = u;
        voxel(otherAxis) = face;
        for (SliceId v = vStart; v < vEnd; ++v) {
          voxel(remainingAxis) = v;
          if (mesh.voxels(voxel.x(), voxel.y(), voxel.z()) != NO_VOXEL) {
            occupiedSlices[u - axisStart][axisListId] = true;
            break;
          }
        }
      }
    }
  };

  processAxis(0);
  processAxis(1);

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

std::vector<RegionId> AxisLayerManager::findCutRegions(RegionId regionId,
                                                       SliceId slice) {
  // std::vector<RegionId> res;
  // bool hasAccessibleParent;
  // do {
  //   hasAccessibleParent = false;
  //   res.push_back(regionId);
  //   Region &region = mesh.regions[regionId];
  //   for (const TraversalCut &cut : regionToChildrenTraversalCuts[regionId])
  //     if (cut.contains(slice)) {
  //       hasAccessibleParent = true;
  //       break;
  //     }
  //   if (hasAccessibleParent) {
  //     RegionLayerId layerId = region.layerIds[axis];
  //     RegionLayerId parent = NO_LAYER;
  //     for (RegionLayerId parents : parentGraph[layerId]) {
  //       if (getAdjacentRegions)
  //     }
  //   }
  // } while (hasAccessibleParent);
  return {};
}
