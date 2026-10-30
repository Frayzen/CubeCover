#include "mesh_loader.hpp"
#include <vtkCell.h>
#include <vtkCellType.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridReader.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <map>
#include <utility>
#include <vector>

bool loadTetrahedralMesh(const std::string &path, TetrahedralMesh &output) {
  constexpr size_t INVALID = std::numeric_limits<size_t>::max();

  // ------------------------------------------------------------
  // 1. Load the legacy VTK file
  // ------------------------------------------------------------

  auto reader = vtkSmartPointer<vtkUnstructuredGridReader>::New();
  reader->SetFileName(path.c_str());
  reader->Update();

  vtkUnstructuredGrid *grid = reader->GetOutput();

  if (!grid)
    return false;

  const vtkIdType nbPoints = grid->GetNumberOfPoints();
  const vtkIdType nbCells = grid->GetNumberOfCells();

  // ------------------------------------------------------------
  // 2. Extract points and tetrahedra only
  // ------------------------------------------------------------

  std::vector<std::array<double, 3>> points(nbPoints);

  for (vtkIdType i = 0; i < nbPoints; ++i) {
    double p[3];
    grid->GetPoint(i, p);
    points[i] = {p[0], p[1], p[2]};
  }

  std::vector<std::array<size_t, 4>> rawTets;

  for (vtkIdType i = 0; i < nbCells; ++i) {
    if (grid->GetCellType(i) != VTK_TETRA)
      continue;

    vtkCell *cell = grid->GetCell(i);

    if (!cell || cell->GetNumberOfPoints() != 4)
      continue;

    std::array<size_t, 4> tet;

    for (size_t j = 0; j < 4; ++j)
      tet[j] = static_cast<size_t>(cell->GetPointId(j));

    rawTets.push_back(tet);
  }

  if (rawTets.empty())
    return false;

  // ------------------------------------------------------------
  // 3. Generate faces and identify boundary faces
  // ------------------------------------------------------------

  struct RawFace {
    std::array<size_t, 3> nodes;
    std::array<ssize_t, 2> adjacentCells = {-1, -1};
  };

  std::vector<RawFace> rawFaces;

  // Canonical face key -> raw face ID
  std::map<std::array<size_t, 3>, size_t> faceMap;

  // Four face IDs for each tetrahedron
  std::vector<std::array<size_t, 4>> rawTetToFaces(rawTets.size());

  // Outward-oriented local faces.
  // The last node of each face is the opposite tetrahedron vertex.
  constexpr size_t localFaces[4][4] = {
      {0, 2, 1, 3}, {0, 1, 3, 2}, {1, 2, 3, 0}, {2, 0, 3, 1}};

  for (size_t t = 0; t < rawTets.size(); ++t) {
    const auto &tet = rawTets[t];

    for (size_t f = 0; f < 4; ++f) {
      size_t a = tet[localFaces[f][0]];
      size_t b = tet[localFaces[f][1]];
      size_t c = tet[localFaces[f][2]];
      size_t d = tet[localFaces[f][3]];

      // Ensure the face normal points away from the
      // opposite vertex of the tetrahedron.

      const auto &A = points[a];
      const auto &B = points[b];
      const auto &C = points[c];
      const auto &D = points[d];

      double ux = B[0] - A[0];
      double uy = B[1] - A[1];
      double uz = B[2] - A[2];

      double vx = C[0] - A[0];
      double vy = C[1] - A[1];
      double vz = C[2] - A[2];

      double nx = uy * vz - uz * vy;
      double ny = uz * vx - ux * vz;
      double nz = ux * vy - uy * vx;

      double dot = nx * (D[0] - A[0]) + ny * (D[1] - A[1]) + nz * (D[2] - A[2]);

      if (dot > 0.0)
        std::swap(b, c);

      std::array<size_t, 3> key = {a, b, c};
      std::sort(key.begin(), key.end());

      auto it = faceMap.find(key);

      size_t faceId;

      if (it == faceMap.end()) {
        faceId = rawFaces.size();

        RawFace face;
        face.nodes = {a, b, c};
        face.adjacentCells = {static_cast<ssize_t>(t), -1};

        rawFaces.push_back(face);
        faceMap[key] = faceId;
      } else {
        faceId = it->second;

        RawFace &face = rawFaces[faceId];

        if (face.adjacentCells[1] != -1)
          return false; // Non-manifold face

        face.adjacentCells[1] = static_cast<ssize_t>(t);
      }

      rawTetToFaces[t][f] = faceId;
    }
  }

  const size_t nbBoundaryFace =
      std::count_if(rawFaces.begin(), rawFaces.end(),
                    [](const RawFace &f) { return f.adjacentCells[1] == -1; });

  // ------------------------------------------------------------
  // 4. Identify boundary vertices and boundary tetrahedra
  // ------------------------------------------------------------

  std::vector<bool> isBoundaryNode(nbPoints, false);
  std::vector<bool> isBoundaryTet(rawTets.size(), false);

  for (const RawFace &face : rawFaces) {
    if (face.adjacentCells[1] != -1)
      continue;

    for (size_t node : face.nodes)
      isBoundaryNode[node] = true;

    isBoundaryTet[face.adjacentCells[0]] = true;
  }

  // Boundary vertices first, then interior vertices.
  std::vector<size_t> nodeOrder;
  nodeOrder.reserve(nbPoints);

  for (size_t i = 0; i < static_cast<size_t>(nbPoints); ++i) {
    if (isBoundaryNode[i])
      nodeOrder.push_back(i);
  }

  const size_t nbBoundaryNode = nodeOrder.size();

  // Only include vertices actually used by retained tetrahedra.
  std::vector<bool> usedNode(nbPoints, false);

  for (const auto &tet : rawTets)
    for (size_t node : tet)
      usedNode[node] = true;

  nodeOrder.clear();

  for (size_t i = 0; i < static_cast<size_t>(nbPoints); ++i) {
    if (usedNode[i] && isBoundaryNode[i])
      nodeOrder.push_back(i);
  }

  const size_t boundaryNodeCount = nodeOrder.size();

  for (size_t i = 0; i < static_cast<size_t>(nbPoints); ++i) {
    if (usedNode[i] && !isBoundaryNode[i])
      nodeOrder.push_back(i);
  }

  std::vector<size_t> nodeRemap(nbPoints, INVALID);

  TetrahedralMesh result;

  result.vertices.reserve(nodeOrder.size());

  for (size_t i = 0; i < nodeOrder.size(); ++i) {
    size_t oldId = nodeOrder[i];
    nodeRemap[oldId] = i;

    const auto &p = points[oldId];
    result.vertices.emplace_back(Vertex{p[0], p[1], p[2]});
  }

  // ------------------------------------------------------------
  // 5. Reorder tetrahedra: boundary first
  // ------------------------------------------------------------

  std::vector<size_t> tetOrder;
  tetOrder.reserve(rawTets.size());

  for (size_t i = 0; i < rawTets.size(); ++i)
    if (isBoundaryTet[i])
      tetOrder.push_back(i);

  result.nbBoundaryTet = tetOrder.size();

  for (size_t i = 0; i < rawTets.size(); ++i)
    if (!isBoundaryTet[i])
      tetOrder.push_back(i);

  std::vector<size_t> tetRemap(rawTets.size());

  result.tets.reserve(rawTets.size());

  for (size_t i = 0; i < tetOrder.size(); ++i) {
    size_t oldId = tetOrder[i];
    tetRemap[oldId] = i;

    std::array<size_t, 4> nodes;

    for (size_t j = 0; j < 4; ++j)
      nodes[j] = nodeRemap[rawTets[oldId][j]];

    result.tets.emplace_back(nodes);
  }

  // ------------------------------------------------------------
  // 6. Reorder faces: boundary first
  // ------------------------------------------------------------

  std::vector<size_t> faceOrder;
  faceOrder.reserve(rawFaces.size());

  for (size_t i = 0; i < rawFaces.size(); ++i)
    if (rawFaces[i].adjacentCells[1] == -1)
      faceOrder.push_back(i);

  result.nbBoundaryFace = faceOrder.size();

  for (size_t i = 0; i < rawFaces.size(); ++i)
    if (rawFaces[i].adjacentCells[1] != -1)
      faceOrder.push_back(i);

  std::vector<size_t> faceRemap(rawFaces.size());

  result.faces.reserve(rawFaces.size());

  for (size_t i = 0; i < faceOrder.size(); ++i) {
    size_t oldId = faceOrder[i];
    faceRemap[oldId] = i;

    const RawFace &face = rawFaces[oldId];

    std::array<size_t, 3> nodes;

    for (size_t j = 0; j < 3; ++j)
      nodes[j] = nodeRemap[face.nodes[j]];

    Face newFace(nodes);

    newFace.adjacentCells[0] = tetRemap[face.adjacentCells[0]];

    newFace.adjacentCells[1] =
        face.adjacentCells[1] == -1
            ? -1
            : static_cast<ssize_t>(tetRemap[face.adjacentCells[1]]);

    result.faces.push_back(newFace);
  }

  // ------------------------------------------------------------
  // 7. Generate edges and identify boundary edges
  // ------------------------------------------------------------

  std::vector<Edge> rawEdges;
  std::map<std::array<size_t, 2>, size_t> edgeMap;

  std::vector<std::array<size_t, 3>> rawFaceToEdges(rawFaces.size());

  for (size_t f = 0; f < rawFaces.size(); ++f) {
    const auto &nodes = rawFaces[f].nodes;

    const size_t localEdges[3][2] = {{0, 1}, {1, 2}, {2, 0}};

    for (size_t e = 0; e < 3; ++e) {
      size_t a = nodes[localEdges[e][0]];
      size_t b = nodes[localEdges[e][1]];

      std::array<size_t, 2> key = {std::min(a, b), std::max(a, b)};

      auto it = edgeMap.find(key);

      size_t edgeId;

      if (it == edgeMap.end()) {
        edgeId = rawEdges.size();

        rawEdges.emplace_back(key);
        edgeMap[key] = edgeId;
      } else {
        edgeId = it->second;
      }

      rawFaceToEdges[f][e] = edgeId;
    }
  }

  // An edge is a boundary edge if it belongs to a boundary face.
  std::vector<bool> isBoundaryEdge(rawEdges.size(), false);

  for (size_t f = 0; f < rawFaces.size(); ++f) {
    if (rawFaces[f].adjacentCells[1] != -1)
      continue;

    for (size_t e = 0; e < 3; ++e)
      isBoundaryEdge[rawFaceToEdges[f][e]] = true;
  }

  // ------------------------------------------------------------
  // 8. Reorder edges: boundary first
  // ------------------------------------------------------------

  std::vector<size_t> edgeOrder;
  edgeOrder.reserve(rawEdges.size());

  for (size_t i = 0; i < rawEdges.size(); ++i)
    if (isBoundaryEdge[i])
      edgeOrder.push_back(i);

  result.nbBoundaryEdge = edgeOrder.size();

  for (size_t i = 0; i < rawEdges.size(); ++i)
    if (!isBoundaryEdge[i])
      edgeOrder.push_back(i);

  std::vector<size_t> edgeRemap(rawEdges.size());

  result.edges.reserve(rawEdges.size());

  for (size_t i = 0; i < edgeOrder.size(); ++i) {
    size_t oldId = edgeOrder[i];
    edgeRemap[oldId] = i;

    std::array<size_t, 2> nodes = {nodeRemap[rawEdges[oldId][0]],
                                   nodeRemap[rawEdges[oldId][1]]};

    result.edges.emplace_back(nodes);
  }

  // ------------------------------------------------------------
  // 9. Build tetToFaces
  // ------------------------------------------------------------

  result.tetToFaces.resize(rawTets.size());

  for (size_t i = 0; i < tetOrder.size(); ++i) {
    size_t oldTetId = tetOrder[i];

    for (size_t f = 0; f < 4; ++f) {
      result.tetToFaces[i][f] =
          static_cast<ssize_t>(faceRemap[rawTetToFaces[oldTetId][f]]);
    }
  }

  // ------------------------------------------------------------
  // 10. Build faceToEdges
  // ------------------------------------------------------------

  result.faceToEdges.resize(rawFaces.size());

  for (size_t i = 0; i < faceOrder.size(); ++i) {
    size_t oldFaceId = faceOrder[i];

    result.faceToEdges[i].resize(3);

    for (size_t e = 0; e < 3; ++e) {
      result.faceToEdges[i][e] = edgeRemap[rawFaceToEdges[oldFaceId][e]];
    }
  }

  // ------------------------------------------------------------
  // 11. Build edgeToBoundaryFaces
  // ------------------------------------------------------------

  result.edgeToBoundaryFaces.resize(rawEdges.size());

  for (auto &adj : result.edgeToBoundaryFaces)
    adj = {-1, -1};

  for (size_t f = 0; f < rawFaces.size(); ++f) {
    if (rawFaces[f].adjacentCells[1] != -1)
      continue;

    size_t newFaceId = faceRemap[f];

    for (size_t e = 0; e < 3; ++e) {
      size_t newEdgeId = edgeRemap[rawFaceToEdges[f][e]];

      auto &adj = result.edgeToBoundaryFaces[newEdgeId];

      if (adj[0] == -1)
        adj[0] = static_cast<ssize_t>(newFaceId);
      else if (adj[1] == -1)
        adj[1] = static_cast<ssize_t>(newFaceId);
    }
  }

  // ------------------------------------------------------------
  // 12. Finalize
  // ------------------------------------------------------------

  result.nbBoundaryNode = boundaryNodeCount;

  output = std::move(result);

  return true;
}
