#ifndef REGION_ITERATOR_HPP
#define REGION_ITERATOR_HPP
#include "ids.hpp"
#include <iterator>

class Region;

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
#endif
