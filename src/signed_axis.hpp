#ifndef SIGNED_AXIS_HPP
#define SIGNED_AXIS_HPP

#include "constants.hpp"
#include "ids.hpp"
#include <cmath>
#include <stdexcept>

using Axis = int;
const Axis NO_AXIS = -1;
enum IntegerSign { POSITIVE, NEGATIVE, ZERO, UNDEFINED };

inline int getIntegerSignValue(const IntegerSign &a) {
  switch (a) {
  case POSITIVE:
    return 1;
  case NEGATIVE:
    return -1;
  case ZERO:
    return 0;
  case UNDEFINED:
    return static_cast<int>(INT_NAN);
  }
  throw std::invalid_argument("Invalid signed axis");
}

inline Vertex getSignedAxisVertex(IntegerSign &sa, Axis a) {
  Vertex v = Vertex::Zero();
  v[a] = getIntegerSignValue(sa);
  return v;
}

using Orientation = char;
enum OrientationBits {
  NO_ORIENTATION = 0,
  TOP = 2 << 1,
  BOTTOM = 2 << 2,
  LEFT = 2 << 3,
  RIGHT = 2 << 4,
  FRONT = 2 << 5,
  BACK = 2 << 6
};

inline OrientationBits orientationBits(Axis a, IntegerSign sign) {
  assert(sign != UNDEFINED);
  if (sign == ZERO)
    return NO_ORIENTATION;
  switch (a) {
  case 0:
    return sign == POSITIVE ? LEFT : RIGHT;
  case 1:
    return sign == POSITIVE ? TOP : BOTTOM;
  case 2:
    return sign == POSITIVE ? FRONT : BACK;
  }
  throw std::runtime_error("Invalid parameters.");
  return NO_ORIENTATION;
}

#endif
