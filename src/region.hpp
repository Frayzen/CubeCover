#ifndef REGION_HPP
#define REGION_HPP

#include "ids.hpp"
#include "region_interface.hpp"
#include "signed_axis.hpp"

#include <array>
#include <ostream>
#include <set>
#include <sys/types.h>

class VoxelMesh;
class RegionInterface;

using SliceMatrix = Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic>;

struct Region {
  void updateBounds();
  std::array<std::vector<RegionInterface>, DIMENSION> children;
  std::array<std::vector<RegionInterface>, DIMENSION> parents;

  using TraversalFn = std::function<bool(Region &, RegionInterface &)>;
  using ProcessingFn = std::function<bool(Region &)>;

  void traverseChildren(const TraversalFn &function,
                        const std::array<bool, DIMENSION> &axisMask,
                        const ProcessingFn &processingFn = nullptr);
  void traverseChildren(const TraversalFn &function, const Axis a) {
    std::array<bool, DIMENSION> axisMask{};
    axisMask[a] = true;
    traverseChildren(function, axisMask);
  };
  void traverseParents(const TraversalFn &function,
                       const std::array<bool, DIMENSION> &axisMask,
                       const ProcessingFn &processingFn = nullptr);
  void traverseParents(const TraversalFn &function, const Axis a) {
    std::array<bool, DIMENSION> axisMask{};
    axisMask[a] = true;
    traverseParents(function, axisMask);
  };

  RegionId id;

  using Voxel = Eigen::Vector3<int>;

  VoxelMesh *meshPtr;

  std::vector<Voxel> voxels;
  Eigen::Vector3<SliceId> startSlices;
  Eigen::Vector3<SliceId> endSlices;

  std::array<IntegerSign, DIMENSION> signedNormals = {ZERO, ZERO, ZERO};

  std::pair<std::vector<Voxel>, std::vector<Voxel>>
  splitByFlagMatrix(const FlagMatrix &flagMatrix, const SliceId sliceId,
                    const Axis axis) const;
  Region(RegionId id, VoxelMesh &mesh) : meshPtr(&mesh), id(id) {}

  auto begin() const { return voxels.begin(); };
  auto end() const { return voxels.end(); };

  bool containsSlice(Axis axis, SliceId sliceId) const;

  void addCutOld(Region &newRegion, Axis axis, SliceId sliceId);
  void addCut(Region &newRegion, Axis axis,
              std::set<RegionId> &targetedNeighbours);

  void updateAdjacency(std::vector<RegionId> oldRegionToNewRegion,
                       const std::array<bool, DIMENSION> &axisMask);
  std::vector<RegionInterface>
  interfaces(const std::array<bool, DIMENSION> &axisMask);

  SliceMatrix getSliceMatrix(Axis a, SliceId sliceId);
  Eigen::Vector2i getCoordsInSlice(Axis a, const Voxel &v) const;
  Voxel getVoxelInSlice(Axis axis, const Eigen::Vector2i &coords,
                        SliceId sliceId) const;
};

inline std::ostream &operator<<(std::ostream &os, const Region &r) {
  os << "REGION [" << r.id << "] : " << std::endl;
  os << "FROM :" << r.startSlices.transpose() << std::endl;
  os << "TO :" << r.endSlices.transpose() << std::endl;
  for (Axis d = 0; d < DIMENSION; d++) {
    os << "CHILDREN [" << d << "] ";
    for (RegionInterface i : r.children[d])
      os << i.regionId << " ";
    os << std::endl;
  }
  for (Axis d = 0; d < DIMENSION; d++) {
    os << "PARENTS [" << d << "] ";
    for (auto i : r.parents[d])
      os << i.regionId << " ";
    os << std::endl;
  }
  return os;
}

#endif
