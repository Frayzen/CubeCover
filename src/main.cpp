#include "constants.hpp"
#include "frame_field.hpp"
#include "mesh_loader.hpp"
#include "mesh_saver.hpp"
#include "process.hpp"
#include "trilaterate.hpp"
#include <Eigen/Core>
#include <Eigen/IterativeLinearSolvers>

#include <Eigen/Eigen>
#include <Eigen/Eigenvalues>
#include <iostream>
#include <unordered_set>

struct BBox {
  Vertex max;
  Vertex min;
};

static size_t getPrincipalAxis(Vertex vector) {
  size_t axisAlignBest;
  double signBest;
  double measureBest = 0;
  for (size_t coord = 0; coord < DIMENSION; coord++) {
    Vertex unit = Vertex::Zero();
    unit[coord] = 1;
    double dotProduct = vector.dot(unit);
    double cur = std::abs(dotProduct);
    const double CONSISTENCY_GAP = 10e-3;
    if (measureBest + CONSISTENCY_GAP < cur) {
      measureBest = cur;
      axisAlignBest = coord;
      signBest = cur;
    }
  }
  return axisAlignBest;
}

static size_t findLocalVertex(const Tetrahedron &tet, size_t globalVertex) {
  for (size_t i = 0; i < NB_NODE_PER_TET; ++i) {
    if (tet[i] == globalVertex)
      return i;
  }

  throw std::runtime_error("Face vertex not found in adjacent tetrahedron");
}

struct DiffTetShapeFunctions {
  Vertex g0;
  Vertex g1;
  Vertex g2;
  Vertex g3;

  Vertex &operator[](size_t i) {
    switch (i) {
    case 0:
      return g0;
    case 1:
      return g1;
    case 2:
      return g2;
    case 3:
      return g3;
    default:
      std::cerr << "Wrong value provided " << i << std::endl;
      assert(false);
    }
  }
};

static double computeVolume(TetrahedralMesh &mesh, size_t tetId) {
  Tetrahedron &currentTet = mesh.tets[tetId];

  Vertex &p0 = mesh.vertices[currentTet[0]];
  Vertex &p1 = mesh.vertices[currentTet[1]];
  Vertex &p2 = mesh.vertices[currentTet[2]];
  Vertex &p3 = mesh.vertices[currentTet[3]];

  Eigen::Vector3d v1(p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]);
  Eigen::Vector3d v2(p2[0] - p0[0], p2[1] - p0[1], p2[2] - p0[2]);
  Eigen::Vector3d v3(p3[0] - p0[0], p3[1] - p0[1], p3[2] - p0[2]);

  Eigen::Matrix3d mat;
  mat.col(0) = v1;
  mat.col(1) = v2;
  mat.col(2) = v3;

  double volume = std::abs(mat.determinant()) / 6.0;

  return volume;
}

static DiffTetShapeFunctions computeDiffShapeFunctions(TetrahedralMesh &mesh,
                                                       size_t tetId) {
  const auto &tet = mesh.tets[tetId];
  const auto &p0 = mesh.vertices[tet[0]];

  Eigen::Matrix3d D;
  for (int j = 1; j < 4; ++j) {
    const auto &p = mesh.vertices[tet[j]];
    for (int k = 0; k < 3; ++k) {
      D(k, j - 1) = p[k] - p0[k];
    }
  }

  const double det = D.determinant();
  if (std::abs(det) < 1e-12) {
    throw std::runtime_error("Degenerate tetrahedron");
  }

  const Eigen::Matrix3d invT = D.inverse().transpose();

  DiffTetShapeFunctions res;
  res.g1 = invT.col(0);
  res.g2 = invT.col(1);
  res.g3 = invT.col(2);
  res.g0 = -(res.g1 + res.g2 + res.g3);

  return res;
}

int main(int argc, char *argv[]) {
  if (argc != 2) {
    std::cerr << "Usage: " << argv[0] << " <mesh-file.vtu|mesh-file.vtk>\n";
    return EXIT_FAILURE;
  }

  TetrahedralMesh mesh;
  if (!loadTetrahedralMesh(argv[1], mesh)) {
    return EXIT_FAILURE;
  }

  size_t nbTets = mesh.tets.size();
  size_t nbConstraints =
      mesh.amountInternalFace * (NB_NODE_PER_FACE - 1) * DIMENSION;

  size_t N = nbTets * NB_UNKNOWN_PER_TET + nbConstraints;
  CubeCoverParameters params(mesh, N);

  std::vector<Eigen::Matrix3d> frames = params.frames;
  fillFrameField(frames, mesh);

  std::vector<Eigen::Vector3d> vertices(mesh.tets.size() * NB_NODE_PER_TET);

  vertices.clear();
  for (size_t i = 0; i < mesh.tets.size(); i++) {
    for (size_t j = 0; j < NB_NODE_PER_TET; j++)
      vertices.push_back(mesh.vertices[mesh.tets[i][j]]);
  }

  saveMeshVtk(vertices, frames, "input.vtk");

  std::cout << "Loaded tetrahedral mesh successfully.\n"
            << "Vertices: " << mesh.vertices.size() << '\n'
            << "Tetrahedra: " << mesh.tets.size() << '\n';

  BBox bbox = {
      .max = {mind, mind, mind},
      .min = {maxd, maxd, maxd},
  };

  for (size_t i = 0; i < mesh.vertices.size(); i++) {
    auto &v = mesh.vertices[i];
    for (size_t coord = 0; coord < DIMENSION; coord++) {
      if (v[coord] < bbox.min[coord])
        bbox.min[coord] = v[coord];
      if (v[coord] > bbox.max[coord])
        bbox.max[coord] = v[coord];
    }
  }

  // The unknowns are in the shape [u1 u2 ... v1 v2 ... w1 w2 ... u5 u6 ... v5
  // v6 ... w5 w6 ... deltas ... ]
  std::cout << "BEGIN SYSTEM BUILD" << std::endl;

  Eigen::SparseMatrix<double> &system = params.system;
  Eigen::VectorX<double> &knowns = params.knowns;
  system.reserve(nbTets * NB_UNKNOWN_PER_TET + nbConstraints);

  // Laplace matrix

  for (size_t i = 0; i < mesh.tets.size(); i++) {
    DiffTetShapeFunctions shapeFuncs = computeDiffShapeFunctions(mesh, i);
    Eigen::Matrix3d frame = frames[i];

    const double vol = computeVolume(mesh, i);

    for (size_t coord = 0; coord < DIMENSION; ++coord) {
      const size_t base = i * NB_UNKNOWN_PER_TET + coord * NB_NODE_PER_TET;

      for (size_t di = 0; di < NB_NODE_PER_TET; ++di) {
        const size_t row = base + di;

        for (size_t dj = 0; dj < NB_NODE_PER_TET; ++dj) {
          const size_t col = base + dj;

          const double coeff = vol * shapeFuncs[di].dot(shapeFuncs[dj]);

          system.coeffRef(row, col) += coeff;
        }

        knowns(row) += vol * shapeFuncs[di].dot(frame.col(coord));
      }
    }
  }

  // Constraints

  auto findCoord = [](Eigen::Matrix3i &pi,
                      size_t coord) -> std::tuple<size_t, int> {
    Eigen::Vector3i direction = Eigen::Vector3i::Zero();
    direction[coord] = 1;
    direction = pi * direction;
    for (size_t i = 0; i < DIMENSION; i++) {
      if (direction[i] != 0)
        return {i, direction[i]};
    }
    std::cerr << "Error in PI " << pi << std::endl;
    assert(false);
  };

  const size_t constraintOffset = nbTets * NB_UNKNOWN_PER_TET;
  size_t constraintId = 0;
  for (size_t faceId = 0; faceId < mesh.faces.size(); faceId++) {
    Face &face = mesh.faces[faceId];
    if (face.adjacentCells[1] == UNDEFINED_ELEMENT)
      continue;

    size_t tetTId = face.adjacentCells[0];
    Tetrahedron &tetT = mesh.tets[tetTId];
    size_t tetSId = face.adjacentCells[1];
    Tetrahedron &tetS = mesh.tets[tetSId];
    Eigen::Matrix3i pi = computePiST(frames[tetTId], frames[tetSId]);

    for (size_t coord = 0; coord < DIMENSION; coord++) {
      size_t colTetT = tetTId * NB_UNKNOWN_PER_TET + coord * NB_NODE_PER_TET;

      auto [coordS, sign] = findCoord(pi, coord);

      size_t colTetS = tetSId * NB_UNKNOWN_PER_TET + coordS * NB_NODE_PER_TET;
      for (size_t nodeId = 0; nodeId < (NB_NODE_PER_FACE - 1); nodeId++) {
        size_t node = face.vertices[nodeId];
        size_t nextId = (nodeId + 1) % NB_NODE_PER_FACE;
        size_t next = face.vertices[nextId];
        size_t row = constraintOffset + constraintId;
        constraintId++;

        // Enforce f(u_t) - f(v_t) - Pi(f(u_s) - f(v_s)) = 0
        knowns(row) = 0;
        // C
        system.insert(row, colTetT + findLocalVertex(tetT, node)) = 1;
        system.insert(row, colTetS + findLocalVertex(tetS, node)) = -sign;
        system.insert(row, colTetT + findLocalVertex(tetT, next)) = -1;
        system.insert(row, colTetS + findLocalVertex(tetS, next)) = sign;
        // C^t
        system.insert(colTetT + findLocalVertex(tetT, node), row) = 1;
        system.insert(colTetS + findLocalVertex(tetS, node), row) = -sign;
        system.insert(colTetT + findLocalVertex(tetT, next), row) = -1;
        system.insert(colTetS + findLocalVertex(tetS, next), row) = sign;
      }
    }
  }
  assert(constraintId == nbConstraints);

  // // Boundary conditions
  // const size_t bcOffset = nbTets * NB_UNKNOWN_PER_TET + nbConstraints;

  // for (size_t coord = 0; coord < DIMENSION; ++coord) {
  //   const size_t row = bcOffset + coord;
  //   const size_t col = coord * NB_NODE_PER_TET; // tet 0, local vertex 0

  //   system.insert(row, col) = 1.0;
  //   system.insert(col, row) = 1.0;
  //   knowns(row) = 0.0;
  // }

  system.makeCompressed();

  std::cout << "END SYSTEM BUILD" << std::endl;

  // SOLVE

  std::vector<size_t> boundaryRows;
  boundaryRows.push_back(0);

  auto getGlobalVertices = [&](size_t tetId, size_t globalNodeId) {
    auto tetNodes = mesh.tets[tetId];
    for (size_t i = 0; i < tetNodes.size(); i++) {
      if (tetNodes[i] == globalNodeId)
        return getLocalVertices(tetId, i);
    }
    assert(false);
  };

  for (size_t faceId = 0; faceId < mesh.faces.size(); faceId++) {
    if (mesh.faces[faceId].adjacentCells[1] != UNDEFINED_ELEMENT)
      continue;
    Face &face = mesh.faces[faceId];
    size_t tetId = face.adjacentCells[0];

    Vertex textureNodes[NB_NODE_PER_FACE];
    for (size_t i = 0; i < NB_NODE_PER_FACE; i++)
      textureNodes[i] = getGlobalVertices(tetId, face.vertices[i]);

    Vertex normal = getNormal(textureNodes);
    size_t axis = getPrincipalAxis(normal);

    std::array<size_t, DIMENSION * NB_NODE_PER_FACE> rows;
    size_t baseTet = NB_UNKNOWN_PER_TET * tetId;
    for (size_t node = 0; node < NB_NODE_PER_FACE; node++) {
      for (size_t d = 0; d < DIMENSION; d++)
        boundaryRows.push_back(baseTet + d * NB_NODE_PER_TET + node);
    }
  }

  size_t nbIter = boundaryRows.size();
  while (!boundaryRows.empty()) {
    std::cout << "---------- ITERATION " << boundaryRows.size() << "/" << nbIter
              << " ----------" << std::endl;
    cubeCoverIteration(params, boundaryRows[boundaryRows.size() - 1]);

    // TODO update here how we set values
    // vertices.clear();
    // auto getLocalVertices = [&](size_t tetId, size_t localNodeId) -> Vertex {
    //   size_t offset = tetId * NB_UNKNOWN_PER_TET;
    //   Vertex res;
    //   for (size_t dim = 0; dim < DIMENSION; dim++)
    //     res[dim] = params.result[offset + dim * NB_NODE_PER_TET +
    //     localNodeId];
    //   return res;
    // };
    // for (size_t i = 0; i < mesh.tets.size(); i++) {
    //   for (size_t j = 0; j < NB_NODE_PER_TET; j++)
    //     vertices.push_back(getLocalVertices(i, j));
    // }
    // saveMeshVtk(vertices, params.frames,
    //             "iteration" + std::to_string(boundaryRows.size()) + ".vtk");

    // boundaryRows.pop_back();
  }
  return EXIT_SUCCESS;
}
