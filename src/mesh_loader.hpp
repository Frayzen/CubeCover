#ifndef MESH_LOADER_HPP
#define MESH_LOADER_HPP

#include "constants.hpp"
#include "vertex.hpp"
#include <cassert>
#include <vtkCell.h>
#include <vtkCellType.h>

#include <array>
#include <string>
#include <vector>

struct Tetrahedron {
  Tetrahedron() {}
  Tetrahedron(std::array<vtkIdType, NB_NODE_PER_TET> v) : v(v) {}
  std::array<vtkIdType, NB_NODE_PER_TET> v;

  constexpr size_t size() const { return NB_NODE_PER_TET; }
  const vtkIdType &operator[](size_t i) const { return v[i]; }
  vtkIdType &operator[](size_t i) { return v[i]; }
  size_t getLocalId(size_t nodeId) const {
    for (size_t i = 0; i < NB_NODE_PER_TET; ++i) {
      if (v[i] == nodeId)
        return i;
    }

    throw std::runtime_error("Face vertex not found in adjacent tetrahedron");
  };

  Vertex getBarycenter(std::vector<Vertex> &meshVertices) const {
    Vertex res = Vertex::Zero();
    for (size_t i = 0; i < NB_NODE_PER_TET; ++i)
      res += meshVertices[v[i]];
    return res / NB_NODE_PER_TET;
  }
};

class TetrahedralMesh;

struct Face {
  std::array<vtkIdType, 3> vertices;
  std::array<ssize_t, 2> adjacentCells;

  Vertex computeNormal(TetrahedralMesh &m);
};

struct TetrahedralMesh {
  std::vector<Vertex> vertices;
  std::vector<Tetrahedron> tets;
  std::vector<Face> faces;
  std::vector<std::array<ssize_t, 4>> tetToFaces;
  size_t amountInternalFace;
};

bool loadTetrahedralMesh(const std::string &path, TetrahedralMesh &output);

#endif
