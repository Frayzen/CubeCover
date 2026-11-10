#ifndef UTILS_HPP
#define UTILS_HPP
#include "ids.hpp"
bool tetrahedronBoxIntersection(const Vertex &A, const Vertex &B,
                                const Vertex &C, const Vertex &D,
                                const Vertex &boxMin, const Vertex &boxMax);
#endif
