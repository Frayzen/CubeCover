#include "utils.hpp"
bool tetrahedronBoxIntersection(const Vertex &A, const Vertex &B,
                                const Vertex &C, const Vertex &D,
                                const Vertex &boxMin, const Vertex &boxMax) {

  // Box center and half-size
  const Vertex boxCenter = (boxMin + boxMax) * 0.5;
  const Vertex boxHalfSize = (boxMax - boxMin) * 0.5;

  // Translate tetrahedron vertices relative to box center
  const Vertex v[4] = {A - boxCenter, B - boxCenter, C - boxCenter,
                       D - boxCenter};

  // Tetrahedron edges
  const int edges[6][2] = {{0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}};

  constexpr double EPSILON = 1e-12;

  // Test a separating axis
  auto separated = [&](const Vertex &axis) -> bool {
    if (axis.squaredNorm() < EPSILON)
      return false;

    // Project tetrahedron onto axis
    double minProj = v[0].dot(axis);
    double maxProj = minProj;

    for (int i = 1; i < 4; ++i) {
      const double p = v[i].dot(axis);
      minProj = std::min(minProj, p);
      maxProj = std::max(maxProj, p);
    }

    // Project box onto axis
    const double radius = boxHalfSize.x() * std::abs(axis.x()) +
                          boxHalfSize.y() * std::abs(axis.y()) +
                          boxHalfSize.z() * std::abs(axis.z());

    // No intersection if projections do not overlap
    return minProj > radius || maxProj < -radius;
  };

  // 1. Box face normals (X, Y, Z)
  if (separated(Vertex::UnitX()) || separated(Vertex::UnitY()) ||
      separated(Vertex::UnitZ()))
    return false;

  // 2. Tetrahedron face normals
  const Vertex faceNormals[4] = {
      (v[1] - v[0]).cross(v[2] - v[0]), (v[1] - v[0]).cross(v[3] - v[0]),
      (v[2] - v[0]).cross(v[3] - v[0]), (v[2] - v[1]).cross(v[3] - v[1])};

  for (const Vertex &axis : faceNormals) {
    if (separated(axis))
      return false;
  }

  // 3. Cross products between tetrahedron edges
  // and box coordinate axes
  const Vertex boxAxes[3] = {Vertex::UnitX(), Vertex::UnitY(), Vertex::UnitZ()};

  for (const auto &edge : edges) {
    const Vertex e = v[edge[1]] - v[edge[0]];

    for (const Vertex &axis : boxAxes) {
      const Vertex crossAxis = e.cross(axis);

      if (separated(crossAxis))
        return false;
    }
  }

  // No separating axis found: intersection exists
  return true;
}
