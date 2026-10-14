#ifndef MESH_LOADER_HPP
#define MESH_LOADER_HPP

#include "vertex.hpp"
#include <cassert>
#include <vtkCell.h>
#include <vtkCellType.h>

#include <array>
#include <string>
#include <vector>

using Tetrahedron = std::array<vtkIdType, 4>;

struct TetrahedralMesh {
  std::vector<Vertex> vertices;
  std::vector<Tetrahedron> tets;
};

bool loadTetrahedralMesh(const std::string &path, TetrahedralMesh &output);

#endif
