#ifndef MESH_SAVER_HPP
#define MESH_SAVER_HPP

#include <Eigen/Dense>
#include <vector>

void saveMesh(const std::vector<Eigen::Vector3d> &vertices,
              const std::string &filename = "output.mesh");
#endif
