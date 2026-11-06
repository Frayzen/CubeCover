#ifndef REGION_HPP
#define REGION_HPP

#include "slice_manager.hpp"
#include "vertex.hpp"
#include <ostream>
#include <sys/types.h>
using RegionId = ssize_t;
const RegionId NO_REGION = -1;

class VoxelMesh;

class Region;
struct Region {
  VoxelMesh *meshPtr;
  Eigen::Vector3<SliceId> startSlices;
  Eigen::Vector3<SliceId> endSlices;

  std::vector<VoxelId> voxels;

  // TODO: might be useful later on
  // std::array<RegionId, DIMENSION> nextRegion;

  Region(VoxelMesh &mesh) : meshPtr(&mesh) {}

  auto begin() const { return voxels.begin(); }

  auto end() const { return voxels.end(); }

  bool containsVoxel(VoxelId voxelId);

  bool containsSlice(Axis axis, SliceId sliceId);
};

inline std::ostream &operator<<(std::ostream &os, const Region &r) {
  os << "REGION : " << std::endl;
  os << "FROM :" << r.startSlices.transpose() << std::endl;
  os << "TO :" << r.endSlices.transpose() << std::endl;
  return os;
}
#endif
