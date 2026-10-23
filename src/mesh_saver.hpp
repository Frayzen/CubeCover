#ifndef MESH_SAVER_HPP
#define MESH_SAVER_HPP

#include <Eigen/Dense>
#include <vector>

void saveMeshGmsh(const std::vector<Eigen::Vector3d> &vertices,
                  const std::string &filename = "output.mesh");
// void saveMeshVtk(const std::vector<Eigen::Vector3d> &vertices,
//                  const std::vector<Eigen::Matrix3d> &frames,
//                  const std::string &path);

void saveMeshCollection(const std::string &path,
                        const std::vector<std::string> &meshFiles);

void saveMeshVtk(const std::vector<Eigen::Vector3d> &vertices,
                 const std::vector<Eigen::Matrix3d> &frames,
                 const std::vector<Eigen::Vector3i> &fixedMasks,
                 const std::string &path);
#endif
