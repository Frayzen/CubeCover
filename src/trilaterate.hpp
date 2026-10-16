#ifndef TRILATERATE_HPP
#define TRILATERATE_HPP

#include <Eigen/Dense>

std::tuple<Eigen::Vector3d, Eigen::Vector3d>
trilaterate3D(const Eigen::Vector3d nodes[3], double distances[3]);
#endif
