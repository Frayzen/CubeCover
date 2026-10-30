#include "mesh_loader.hpp"
#include "constants.hpp"

#include <Eigen/Eigen>
#include <array>
#include <stdexcept>
#include <vector>

#include <vtkCell.h>
#include <vtkCellType.h>
#include <vtkErrorCode.h>
#include <vtkPoints.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridReader.h>
#include <vtkXMLUnstructuredGridReader.h>

// ------------------------------------------------------------
// Face normal
// ------------------------------------------------------------

Vertex Face::computeNormal(TetrahedralMesh &m) {
  const Vertex &p0 = m.vertices[nodes[0]];
  const Vertex &p1 = m.vertices[nodes[1]];
  const Vertex &p2 = m.vertices[nodes[2]];

  Vertex normal = (p1 - p0).cross(p2 - p0);

  const double norm = normal.norm();

  if (norm < 1e-12)
    throw std::runtime_error("Degenerate face: cannot compute normal");

  normal /= norm;

  const Vertex faceCenter = (p0 + p1 + p2) / 3.0;

  const ssize_t cellId = adjacentCells[0];

  if (cellId == UNDEFINED_ELEMENT)
    throw std::runtime_error("Face has no adjacent tetrahedron");

  const Vertex tetCenter = m.tets[cellId].getBarycenter(m.vertices);

  if (normal.dot(tetCenter - faceCenter) > 0.0)
    normal = -normal;

  return normal;
}

// ------------------------------------------------------------
// Tetrahedron barycenter
// ------------------------------------------------------------

Vertex Tetrahedron::getBarycenter(std::vector<Vertex> &meshVertices) const {
  Vertex res = Vertex::Zero();

  for (size_t i = 0; i < NB_NODE_PER_TET; ++i)
    res += meshVertices[nodes[i]];

  return res / static_cast<double>(NB_NODE_PER_TET);
}

// ------------------------------------------------------------
// Edge center
// ------------------------------------------------------------

Vertex Edge::computeCenter(TetrahedralMesh &m) const {
  Vertex res = Vertex::Zero();

  for (size_t i = 0; i < nodes.size(); ++i)
    res += m.vertices[nodes[i]];

  return res / static_cast<double>(nodes.size());
}

// ------------------------------------------------------------
// Edge direction
// ------------------------------------------------------------

Vertex Edge::computeDir(TetrahedralMesh &m) const {
  return m.vertices[nodes[0]] - m.vertices[nodes[1]];
}

bool Edge::containsNode(size_t nodeId) const {
  return (nodes[0] == nodeId) || (nodes[1] == nodeId);
};
