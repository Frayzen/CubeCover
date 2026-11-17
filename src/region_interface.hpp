#ifndef REGION_INTERFACE_HPP
#define REGION_INTERFACE_HPP

class Region;

#include "old/signed_axis.hpp"
struct RegionInterface {

  IntegerSign sign;
  Axis normal;
  Axis nextAxis;
  Axis prevAxis;

  struct AllowedCuts {
    SliceId from;
    SliceId to;
    AllowedCuts(SliceId from, SliceId to) : from(from), to(to) {};
    inline bool contains(SliceId val) const { return val >= from && val <= to; }
  };

  std::vector<AllowedCuts> rangesNextAxis;
  std::vector<AllowedCuts> rangesPrevAxis;

  std::vector<AllowedCuts> &get(Axis a) {
    assert(a != normal);
    if (a == nextAxis)
      return rangesNextAxis;
    assert(a == prevAxis);
    return rangesPrevAxis;
  }

  const std::vector<AllowedCuts> &get(Axis a) const {
    assert(a != normal);
    if (a == nextAxis)
      return rangesNextAxis;
    assert(a == prevAxis);
    return rangesPrevAxis;
  }

  inline bool contains(Axis a, SliceId val) const {
    for (AllowedCuts r : get(a))
      if (r.contains(val))
        return true;
    return false;
  }
  RegionInterface(Axis normal, IntegerSign sign, RegionId regionId)
      : regionId(regionId), normal(normal), sign(sign),
        nextAxis((normal + 1) % DIMENSION),
        prevAxis((normal + 2) % DIMENSION) {};
  RegionId regionId;

  void updateRanges(Region &first, Region &second);
};
#endif
