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
struct Face {
  std::array<vtkIdType, 3> vertices;
  std::array<vtkIdType, 2> adjacentCells;
};

struct TetrahedralMesh {
  std::vector<Vertex> vertices;
  std::vector<Tetrahedron> tets;
  std::vector<Face> faces;
  size_t amountInternalFace;
};

bool loadTetrahedralMesh(const std::string &path, TetrahedralMesh &output);

#endif
