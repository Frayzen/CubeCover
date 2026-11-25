#ifndef REGION_INTERFACE_HPP
#define REGION_INTERFACE_HPP

class DualRegion;

enum InterfacePosition { IP_TOP, IP_BOT, IP_LEFT, IP_RIGHT };

#include "signed_axis.hpp"

struct DualRegionInterface {

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

  DualRegionInterface(Axis normal, IntegerSign sign, DualRegionId regionId)
      : regionId(regionId), normal(normal), sign(sign),
        nextAxis((normal + 1) % DIMENSION),
        prevAxis((normal + 2) % DIMENSION) {};
  DualRegionId regionId;

  InterfacePosition getInterfacePosition(Axis cutDirection) const;
};

bool findPathInMatrix(
    const Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> &mat,
    const InterfacePosition &from, const InterfacePosition &to);

bool findPathInMatrixFrom(
    const Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> &mat,
    const Eigen::Vector2i &from, const InterfacePosition &to);
#endif
