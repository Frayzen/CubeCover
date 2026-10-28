#ifndef MESH_SAVER_HPP
#define MESH_SAVER_HPP

#include "mesh_loader.hpp"
#include "vertex.hpp"
#include <Eigen/Dense>
#include <vector>
#include <vtkCell.h>

void saveMeshGmsh(const std::vector<Eigen::Vector3d> &vertices,
                  const std::string &filename = "output.mesh");
// void saveMeshVtk(const std::vector<Eigen::Vector3d> &vertices,
//                  const std::vector<Eigen::Matrix3d> &frames,
//                  const std::string &path);

void saveMeshVtk(const std::vector<Vertex> &vertices,
                 const std::vector<Tetrahedron> &tets,
                 const std::vector<Eigen::Matrix3d> &frames,
                 const std::string &path);
#endif
