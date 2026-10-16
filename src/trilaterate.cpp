#include "trilaterate.hpp"
#include <cmath>

std::tuple<Eigen::Vector3d, Eigen::Vector3d>
trilaterate3D(const Eigen::Vector3d nodes[3], double distances[3]) {
  auto &p1 = nodes[0];
  auto &p2 = nodes[1];
  auto &p3 = nodes[2];

  auto &r1 = distances[0];
  auto &r2 = distances[1];
  auto &r3 = distances[2];

  Eigen::Vector3d p21 = p2 - p1;
  double d = p21.norm();
  if (d < 1e-9) {
    throw std::runtime_error("p1 and p2 cannot be the same point.");
  }
  Eigen::Vector3d ex = p21 / d;

  Eigen::Vector3d p31 = p3 - p1;
  double i = ex.dot(p31);

  Eigen::Vector3d aux = p31 - i * ex;
  double j = aux.norm();
  if (j < 1e-9) {
    throw std::runtime_error("p1, p2, and p3 are collinear.");
  }
  Eigen::Vector3d ey = aux / j;
  Eigen::Vector3d ez = ex.cross(ey);

  // 2. Solve for x and y in the local coordinate system
  double x = (r1 * r1 - r2 * r2 + d * d) / (2.0 * d);
  double y = (r1 * r1 - r3 * r3 + i * i + j * j - 2.0 * i * x) / (2.0 * j);

  // 3. Solve for z^2
  double z_squared = r1 * r1 - x * x - y * y;

  if (z_squared < -1e-9) {
    throw std::runtime_error("No real intersection");
  }

  double z = (z_squared > 0.0) ? std::sqrt(z_squared) : 0.0;

  Eigen::Vector3d sol1 = p1 + x * ex + y * ey + z * ez;
  Eigen::Vector3d sol2 = p1 + x * ex + y * ey - z * ez;
  return {sol1, sol2};
}
