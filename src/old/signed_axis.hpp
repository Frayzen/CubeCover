#ifndef SIGNED_AXIS_HPP
#define SIGNED_AXIS_HPP

#include "constants.hpp"
#include "vertex.hpp"
#include <cmath>
#include <stdexcept>

using Axis = int;
enum IntegerSign { POSITIVE, NEGATIVE, ZERO, UNDEFINED };

inline int getIntegerSignValue(IntegerSign &a) {
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

#endif
