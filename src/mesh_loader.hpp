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

class TetrahedralMesh;

struct Edge {
  std::array<size_t, 2> nodes;
  Edge() {};
  Edge(const std::array<size_t, 2> &nodes) : nodes(nodes) {}
  const size_t &operator[](size_t i) const { return nodes[i]; }
  size_t &operator[](size_t i) { return nodes[i]; }
  Vertex computeCenter(TetrahedralMesh &m) const;
  Vertex computeDir(TetrahedralMesh &m) const;
  bool containsNode(size_t nodeId) const;
};

struct Tetrahedron {
  Tetrahedron() {}
  Tetrahedron(const std::array<size_t, NB_NODE_PER_TET> v) : nodes(v) {}
  std::array<size_t, NB_NODE_PER_TET> nodes;

  const size_t &operator[](size_t i) const { return nodes[i]; }
  size_t &operator[](size_t i) { return nodes[i]; }
  size_t getLocalId(size_t nodeId) const {
    for (size_t i = 0; i < NB_NODE_PER_TET; ++i) {
      if (nodes[i] == nodeId)
        return i;
    }
    throw std::runtime_error("Face vertex not found in adjacent tetrahedron");
  };

  Vertex getBarycenter(std::vector<Vertex> &meshVertices) const;
};

struct Face {
  std::array<size_t, 3> nodes;
  Face(const std::array<size_t, 3> nodes) : nodes(nodes) {}
  std::array<ssize_t, 2> adjacentCells;

  Vertex computeNormal(TetrahedralMesh &m);
};

struct TetrahedralMesh {
  std::vector<Vertex> vertices;  // boundary vertices at the start
  std::vector<Tetrahedron> tets; // boundary tets at the start
  std::vector<Face> faces;       // boudnary faces at the start
  std::vector<Edge> edges;       // boundary edges at the start
  size_t nbBoundaryTet;
  size_t nbBoundaryFace;
  size_t nbBoundaryNode;
  size_t nbBoundaryEdge;
  std::vector<std::array<ssize_t, 4>> tetToFaces;
  std::vector<std::vector<size_t>> faceToEdges;
  std::vector<std::array<ssize_t, 2>> edgeToBoundaryFaces;
  std::vector<Vertex> getBoundaryVertices() {
    return std::vector<Vertex>(vertices.begin(),
                               vertices.begin() + nbBoundaryNode);
  }
};

bool loadTetrahedralMesh(const std::string &path, TetrahedralMesh &output);

#endif
