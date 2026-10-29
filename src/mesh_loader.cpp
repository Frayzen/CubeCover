#include "mesh_loader.hpp"
#include "constants.hpp"

#include <Eigen/Eigen>
#include <algorithm>
#include <array>
#include <filesystem>
#include <iostream>
#include <map>
#include <stdexcept>
#include <utility>

#include <vtkCell.h>
#include <vtkCellType.h>
#include <vtkErrorCode.h>
#include <vtkPoints.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridReader.h>
#include <vtkXMLUnstructuredGridReader.h>

// ------------------------------------------------------------
// Read VTK grid
// ------------------------------------------------------------

static vtkSmartPointer<vtkUnstructuredGrid> readGrid(const std::string &path) {
  const std::filesystem::path filePath(path);
  const std::string extension = filePath.extension().string();

  if (extension == ".vtu" || extension == ".VTU") {
    auto reader = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();

    reader->SetFileName(path.c_str());
    reader->Update();

    if (reader->GetErrorCode() != vtkErrorCode::NoError)
      return nullptr;

    auto grid = vtkSmartPointer<vtkUnstructuredGrid>::New();
    grid->DeepCopy(reader->GetOutput());

    return grid;
  }

  if (extension == ".vtk" || extension == ".VTK") {
    auto reader = vtkSmartPointer<vtkUnstructuredGridReader>::New();

    reader->SetFileName(path.c_str());
    reader->Update();

    if (reader->GetErrorCode() != vtkErrorCode::NoError)
      return nullptr;

    auto grid = vtkSmartPointer<vtkUnstructuredGrid>::New();
    grid->DeepCopy(reader->GetOutput());

    return grid;
  }

  std::cerr << "Unsupported mesh extension: " << extension << '\n';
  return nullptr;
}

// ------------------------------------------------------------
// Load tetrahedral mesh
// ------------------------------------------------------------

static void processEdges(TetrahedralMesh &candidate, bool boundaryEdges) {
  std::map<std::array<vtkIdType, 2>, size_t> edgeMap;
  for (size_t faceId = 0; faceId < candidate.faces.size(); ++faceId) {
    bool isBoundary =
        candidate.faces[faceId].adjacentCells[1] == UNDEFINED_ELEMENT;
    if (boundaryEdges && !isBoundary)
      continue;
    if (!boundaryEdges && isBoundary)
      continue;

    const Face &face = candidate.faces[faceId];

    auto &faceEdges = candidate.faceToEdges[faceId];
    faceEdges.reserve(3);

    // Three edges per triangular face.
    constexpr std::array<std::array<size_t, 2>, 3> localEdges = {
        {{{0, 1}}, {{1, 2}}, {{2, 0}}}};

    for (const auto &localEdge : localEdges) {
      vtkIdType a = face.vertices[localEdge[0]];
      vtkIdType b = face.vertices[localEdge[1]];

      std::array<vtkIdType, 2> key = {std::min(a, b), std::max(a, b)};

      auto it = edgeMap.find(key);

      size_t edgeId;

      if (it == edgeMap.end()) {
        edgeId = candidate.edges.size();

        Edge edge;
        edge.v = key;

        candidate.edges.push_back(edge);
        edgeMap.emplace(key, edgeId);
      } else {
        edgeId = it->second;
      }

      faceEdges.push_back(edgeId);
    }
  }
}

bool loadTetrahedralMesh(const std::string &path, TetrahedralMesh &output) {
  auto grid = readGrid(path);

  if (!grid) {
    std::cerr << "Failed to read mesh: " << path << '\n';
    return false;
  }

  const vtkIdType nbPoints = grid->GetNumberOfPoints();
  const vtkIdType nbCells = grid->GetNumberOfCells();

  if (nbPoints == 0 || nbCells == 0) {
    std::cerr << "Mesh contains no points or cells.\n";
    return false;
  }

  // --------------------------------------------------------
  // Validate input
  // --------------------------------------------------------

  for (vtkIdType cellId = 0; cellId < nbCells; ++cellId) {
    if (grid->GetCellType(cellId) != VTK_TETRA) {
      std::cerr << "Cell " << cellId << " is not a tetrahedron.\n";
      return false;
    }

    vtkCell *cell = grid->GetCell(cellId);

    if (!cell || cell->GetNumberOfPoints() != NB_NODE_PER_TET ||
        cell->GetNumberOfFaces() != 4) {
      std::cerr << "Invalid tetrahedron at cell " << cellId << ".\n";
      return false;
    }
  }

  TetrahedralMesh candidate;
  candidate.amountInternalFace = 0;

  // --------------------------------------------------------
  // Load vertices
  // --------------------------------------------------------

  candidate.vertices.reserve(static_cast<size_t>(nbPoints));

  vtkPoints *points = grid->GetPoints();

  for (vtkIdType pointId = 0; pointId < nbPoints; ++pointId) {
    double p[3];
    points->GetPoint(pointId, p);

    candidate.vertices.emplace_back(Vertex{p[0], p[1], p[2]});
  }

  // --------------------------------------------------------
  // Load tetrahedra
  // --------------------------------------------------------

  candidate.tets.reserve(static_cast<size_t>(nbCells));
  candidate.tetToFaces.resize(static_cast<size_t>(nbCells));

  // Initialize tet-to-face connectivity.
  for (auto &faces : candidate.tetToFaces) {
    faces.fill(UNDEFINED_ELEMENT);
  }

  // Maps a sorted face vertex triplet to its global face ID.
  std::map<std::array<vtkIdType, 3>, size_t> faceMap;

  for (vtkIdType cellId = 0; cellId < nbCells; ++cellId) {
    vtkCell *cell = grid->GetCell(cellId);

    Tetrahedron tet;

    for (size_t i = 0; i < NB_NODE_PER_TET; ++i) {
      const vtkIdType pointId = cell->GetPointId(i);

      if (pointId < 0 || pointId >= nbPoints) {
        std::cerr << "Invalid vertex index in tetrahedron " << cellId << ".\n";
        return false;
      }

      tet[i] = pointId;
    }

    candidate.tets.push_back(tet);

    // ----------------------------------------------------
    // Extract the four faces of this tetrahedron.
    // ----------------------------------------------------

    for (size_t localFaceId = 0; localFaceId < 4; ++localFaceId) {
      vtkCell *vtkFace = cell->GetFace(localFaceId);

      if (!vtkFace || vtkFace->GetNumberOfPoints() != 3) {
        std::cerr << "Invalid face in tetrahedron " << cellId << ".\n";
        return false;
      }

      std::array<vtkIdType, 3> key = {vtkFace->GetPointId(0),
                                      vtkFace->GetPointId(1),
                                      vtkFace->GetPointId(2)};

      // The key is independent of face orientation.
      std::sort(key.begin(), key.end());

      auto it = faceMap.find(key);

      size_t faceId;

      if (it == faceMap.end()) {
        // First tetrahedron adjacent to this face.
        faceId = candidate.faces.size();

        Face face;
        face.vertices = key;
        face.adjacentCells = {static_cast<ssize_t>(cellId), UNDEFINED_ELEMENT};

        candidate.faces.push_back(face);
        faceMap.emplace(key, faceId);
      } else {
        faceId = it->second;

        Face &face = candidate.faces[faceId];

        if (face.adjacentCells[1] != UNDEFINED_ELEMENT) {
          std::cerr << "Non-manifold face shared by more "
                       "than two tetrahedra.\n";
          return false;
        }

        // Second tetrahedron adjacent to this face.
        face.adjacentCells[1] = static_cast<ssize_t>(cellId);

        candidate.amountInternalFace++;
      }

      // Preserve the local face ordering of the tetrahedron.
      candidate.tetToFaces[cellId][localFaceId] = static_cast<ssize_t>(faceId);
    }
  }

  // --------------------------------------------------------
  // Build unique edges and face-to-edge connectivity
  // --------------------------------------------------------

  // Maps a sorted vertex pair to its global edge ID.

  candidate.faceToEdges.resize(candidate.faces.size());

  processEdges(candidate, true);
  candidate.nbBoundaryEdge = candidate.edges.size();
  processEdges(candidate, false);

  // --------------------------------------------------------
  // Commit result
  // --------------------------------------------------------

  output = std::move(candidate);

  return true;
}

// ------------------------------------------------------------
// Face normal
// ------------------------------------------------------------

Vertex Face::computeNormal(TetrahedralMesh &m) {
  const Vertex &p0 = m.vertices[vertices[0]];
  const Vertex &p1 = m.vertices[vertices[1]];
  const Vertex &p2 = m.vertices[vertices[2]];

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
    res += meshVertices[v[i]];

  return res / static_cast<double>(NB_NODE_PER_TET);
}

Vertex Edge::computeCenter(TetrahedralMesh &m) const {
  Vertex res = Vertex::Zero();
  for (size_t i = 0; i < v.size(); ++i)
    res += m.vertices[v[i]];
  return res / static_cast<double>(v.size());
};

Vertex Edge::computeDir(TetrahedralMesh &m) const {
  return m.vertices[v[0]] - m.vertices[v[1]];
};
