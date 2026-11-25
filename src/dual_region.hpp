#ifndef DUAL_REGION_HPP
#define DUAL_REGION_HPP

#include "dual_region_interface.hpp"
#include "ids.hpp"
#include "signed_axis.hpp"

#include <array>
#include <ostream>
#include <sys/types.h>

class VoxelMesh;
class DualRegionInterface;

using SliceMatrix = Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic>;

struct DualRegion {
  void updateBounds();
  std::array<std::vector<DualRegionInterface>, DIMENSION> children;
  // std::array<std::vector<RegionInterface>, DIMENSION> parents;

  using TraversalFn = std::function<bool(DualRegion &, DualRegionInterface &)>;
  using ProcessingFn = std::function<bool(DualRegion &)>;

  DualRegionId id;

  using Voxel = Eigen::Vector3<int>;

  VoxelMesh *meshPtr;

  std::vector<Voxel> voxels;
  Eigen::Vector3<SliceId> startSlices;
  Eigen::Vector3<SliceId> endSlices;

  // those normals is ZERO everywhere but on the extermity of the region graph
  std::array<IntegerSign, DIMENSION> boundaryNormals = {ZERO, ZERO, ZERO};
  // those normals is never zero and represents the existing faces of the reigon
  std::array<IntegerSign, DIMENSION> facesNormals = {ZERO, ZERO, ZERO};

  std::array<SurfaceId, DIMENSION> surfaceIds;

  std::pair<std::vector<Voxel>, std::vector<Voxel>>
  splitByFlagMatrix(const FlagMatrix &flagMatrix, const SliceId sliceId,
                    const Axis axis) const;
  DualRegion(DualRegionId id, VoxelMesh &mesh) : meshPtr(&mesh), id(id) {}

  auto begin() const { return voxels.begin(); };
  auto end() const { return voxels.end(); };

  bool containsSlice(Axis axis, SliceId sliceId) const;

  std::vector<DualRegionInterface *>
  interfaces(const std::array<bool, DIMENSION> &axisMask = {true, true, true});

  std::vector<DualRegionInterface *> interfaces(Axis axis) {
    std::array<bool, DIMENSION> axisMask = {};
    axisMask[axis] = true;
    return interfaces(axisMask);
  }

  SliceMatrix getSliceMatrix(Axis a, SliceId sliceId);
  Eigen::Vector2i getCoordsInSlice(Axis a, const Voxel &v) const;
  Voxel getVoxelInSlice(Axis axis, const Eigen::Vector2i &coords,
                        SliceId sliceId) const;
};

inline std::ostream &operator<<(std::ostream &os, const DualRegion &r) {
  os << "REGION [" << r.id << "] : " << std::endl;
  os << "FROM :" << r.startSlices.transpose() << std::endl;
  os << "TO :" << r.endSlices.transpose() << std::endl;
  for (Axis d = 0; d < DIMENSION; d++) {
    os << "CHILDREN [" << d << "] ";
    for (DualRegionInterface i : r.children[d])
      os << i.regionId << " ";
    os << std::endl;
  }
  // for (Axis d = 0; d < DIMENSION; d++) {
  // os << "PARENTS [" << d << "] ";
  // for (auto i : r.parents[d])
  //   os << i.regionId << " ";
  // os << std::endl;
  // }
  return os;
}

#endif
