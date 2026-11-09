#ifndef REGION_HPP
#define REGION_HPP

#include "signed_axis.hpp"
#include "slice_manager.hpp"
#include "vertex.hpp"
#include <ostream>
#include <sys/types.h>
using RegionId = ssize_t;
const RegionId NO_REGION = -1;

class VoxelMesh;

using RegionLayerId = size_t;

class Region;
struct Region {
  VoxelMesh *meshPtr;
  Eigen::Vector3<SliceId> startSlices;
  Eigen::Vector3<SliceId> endSlices;

  std::array<IntegerSign, DIMENSION> signedNormals = {ZERO, ZERO, ZERO};

  std::vector<VoxelId> voxels;

  std::array<RegionLayerId, DIMENSION> layerIds;

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
