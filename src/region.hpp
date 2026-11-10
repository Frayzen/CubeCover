#ifndef REGION_HPP
#define REGION_HPP

#include "ids.hpp"
#include "signed_axis.hpp"

#include <array>
#include <cstddef>
#include <iterator>
#include <ostream>
#include <sys/types.h>

class VoxelMesh;

struct Region {
  using Voxel = Eigen::Vector3<int>;

  class RegionIterator {
  public:
    using iterator_category = std::forward_iterator_tag;
    using value_type = Voxel;
    using difference_type = std::ptrdiff_t;
    using pointer = void;
    using reference = Voxel;

    RegionIterator(const Region *region, bool atEnd = false);

    Voxel operator*() const;

    RegionIterator &operator++();

    RegionIterator operator++(int) {
      RegionIterator copy = *this;
      ++(*this);
      return copy;
    }

    bool operator==(const RegionIterator &other) const;
    bool operator!=(const RegionIterator &other) const {
      return !(*this == other);
    }

  private:
    const Region *region_;
    int x_ = 0;
    int y_ = 0;
    int z_ = 0;
    bool atEnd_ = true;

    void advanceToOccupied();
  };

  VoxelMesh *meshPtr;
  Eigen::Vector3<SliceId> startSlices;
  Eigen::Vector3<SliceId> endSlices;

  std::array<IntegerSign, DIMENSION> signedNormals = {ZERO, ZERO, ZERO};
  std::array<RegionLayerId, DIMENSION> layerIds;

  Region(VoxelMesh &mesh) : meshPtr(&mesh) {}

  RegionIterator begin() const { return RegionIterator(this); }

  RegionIterator end() const { return RegionIterator(this, true); }

  bool containsSlice(Axis axis, SliceId sliceId);
};

inline std::ostream &operator<<(std::ostream &os, const Region &r) {
  os << "REGION : " << std::endl;
  os << "FROM :" << r.startSlices.transpose() << std::endl;
  os << "TO :" << r.endSlices.transpose() << std::endl;
  return os;
}

#endif
