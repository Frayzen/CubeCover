#ifndef AXIS_LAYER_MANAGER_HPP
#define AXIS_LAYER_MANAGER_HPP

#include "graph.hpp"
#include "ids.hpp"
#include "region.hpp"

class VoxelMesh;

// a TraversalCut describes from what layer to what layer of the face the cut is
// "passing". In other worlds, it describes the slice that connect two regions
struct TraversalCut {
  SliceId from;
  SliceId to;
  inline bool contains(SliceId val) { return val >= from && val <= to; }
};

using TraversalCuts = std::vector<TraversalCut>;

class AxisLayerManager {
private:
  Axis axis;

  // A node here is a RegionLayerId
  Graph childGraph;
  Graph parentGraph;
  VoxelMesh &mesh;

  std::vector<TraversalCuts> regionToChildrenTraversalCuts;
  std::vector<std::vector<RegionId>> layerToRegions;

  void updateTraversalCuts(RegionId regionId);

public:
  AxisLayerManager(VoxelMesh &mesh, Axis axis)
      : mesh(mesh), axis(axis), childGraph(0), parentGraph(0) {};

  Graph &graph() { return childGraph; };
  RegionLayerId prepareSplit(RegionLayerId layerId);
  void concretizeSplit(RegionLayerId layerId, std::vector<RegionId> &regions);
  void initialize(Region &initRegion);
  void splitRegionFromAnotherAxis(RegionId splittedRegionId, RegionId regionId);

  bool isCutAllowed(RegionLayerId id, SliceId slice) const;

  friend std::ostream &operator<<(std::ostream &os, const AxisLayerManager &m);
};

inline std::ostream &operator<<(std::ostream &os, const AxisLayerManager &m) {
  os << "LAYER MANAGER " << m.axis << " : " << std::endl;
  RegionId rId = 0;
  for (auto cuts : m.regionToChildrenTraversalCuts) {
    os << "[REGION " << rId << "]" << std::endl;
    for (auto cut : cuts) {
      os << cut.from << " - " << cut.to << std::endl;
    }
    rId++;
  }

  return os;
}
#endif
