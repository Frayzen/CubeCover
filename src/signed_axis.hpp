#ifndef SIGNED_AXIS_HPP
#define SIGNED_AXIS_HPP

#include "vertex.hpp"
#include <cmath>
#include <cstddef>
#include <stdexcept>

using Axis = size_t;
enum SignedAxis { POSITIVE, NEGATIVE, ZERO, UNDEFINED };

inline int getSignedAxisValue(SignedAxis &a) {
  switch (a) {
  case POSITIVE:
    return 1;
  case NEGATIVE:
    return -1;
  case ZERO:
    return 0;
  case UNDEFINED:
    return static_cast<int>(NAN);
  }
  throw std::invalid_argument("Invalid signed axis");
}

inline Vertex getSignedAxisVertex(SignedAxis &sa, Axis a) {
  Vertex v = Vertex::Zero();
  v[a] = getSignedAxisValue(sa);
  return v;
}

#endif
