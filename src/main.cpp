#include "constants.hpp"
#include "mesh_loader.hpp"
#include "mesh_saver.hpp"
#include "span_tree.hpp"
#include "trilaterate.hpp"
#include <Eigen/Core>
// #include <Eigen/SparseLU>
#include <Eigen/IterativeLinearSolvers>

#include <Eigen/Eigen>
#include <Eigen/Eigenvalues>
#include <iostream>

struct BBox {
  Vertex max;
  Vertex min;
};

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
  std::vector<size_t> spanTree = convertToSpanTree(mesh);

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

  size_t nbTets = mesh.tets.size();
  size_t nbConstraints = mesh.amountInternalFace * NB_NODE_PER_FACE * DIMENSION;

  size_t bcAmount = DIMENSION * 2;
  size_t N = nbTets * NB_UNKNOWN_PER_TET + nbConstraints + bcAmount;

  // The unknowns are in the shape [u1 u2 ... v1 v2 ... w1 w2 ... u5 u6 ... v5
  // v6 ... w5 w6 ... deltas ... ]
  Eigen::SparseMatrix<double> system(N, N);
  system.reserve(nbTets * NB_UNKNOWN_PER_TET + nbConstraints);

  Eigen::VectorX<double> knowns = Eigen::VectorXd::Zero(N);

  // Laplace matrix

  for (size_t i = 0; i < mesh.tets.size(); i++) {
    DiffTetShapeFunctions shapeFuncs = computeDiffShapeFunctions(mesh, i);
    Eigen::Matrix3d frame;   // TODO set the frame not constant
    frame.col(0) << 1, 0, 0; // U
    frame.col(1) << 0, 1, 0; // V
    frame.col(2) << 0, 0, 1; // W

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
    // Here, we assume that Pi is Identity
    for (size_t coord = 0; coord < DIMENSION; coord++) {
      size_t colTetT = tetTId * NB_UNKNOWN_PER_TET + coord * NB_NODE_PER_TET;
      size_t colTetS = tetSId * NB_UNKNOWN_PER_TET + coord * NB_NODE_PER_TET;
      for (size_t nodeId = 0; nodeId < NB_NODE_PER_FACE; nodeId++) {
        size_t node = face.vertices[nodeId];
        size_t nextId = (nodeId + 1) % NB_NODE_PER_FACE;
        size_t next = face.vertices[nextId];
        size_t row = constraintOffset + constraintId;
        constraintId++;

        // Enforce f(u_t) - f(v_t) - Pi(f(u_s) - f(v_s)) = 0
        knowns(row) = 0;
        // C
        system.insert(row, colTetT + findLocalVertex(tetT, node)) = 1;
        system.insert(row, colTetS + findLocalVertex(tetS, node)) =
            -1; // TODO apply Pi
        system.insert(row, colTetT + findLocalVertex(tetT, next)) = -1;
        system.insert(row, colTetS + findLocalVertex(tetS, next)) =
            1; // TODO apply Pi
        // C^t
        system.insert(colTetT + findLocalVertex(tetT, node), row) = 1;
        system.insert(colTetS + findLocalVertex(tetS, node), row) =
            -1; // TODO apply Pi
        system.insert(colTetT + findLocalVertex(tetT, next), row) = -1;
        system.insert(colTetS + findLocalVertex(tetS, next), row) =
            1; // TODO apply Pi
      }

      // // Enforce f(u_t)
      // size_t row = constraintOffset + constraintId;
      // constraintId++;
      // knowns(row) = 0;
      // // C
      // system.insert(row, colTetT) = 1;
      // system.insert(row, colTetS) = -1; // TODO apply Pi
      // // C^t
      // system.insert(colTetT, row) = 1;
      // system.insert(colTetS, row) = -1; // TODO apply Pi
    }
  }
  assert(constraintId == nbConstraints);

  // Boundary conditions
  for (size_t coord = 0; coord < DIMENSION; coord++) {
    system.insert(nbTets * NB_UNKNOWN_PER_TET + nbConstraints + coord,
                  coord * NB_NODE_PER_TET) = 1;
    system.insert(nbTets * NB_UNKNOWN_PER_TET + nbConstraints + coord,
                  coord * NB_NODE_PER_TET + 1) = 1;
    // C^t
    system.insert(coord * NB_NODE_PER_TET,
                  nbTets * NB_UNKNOWN_PER_TET + nbConstraints + coord) = 1;
    system.insert(coord * NB_NODE_PER_TET + 1,
                  nbTets * NB_UNKNOWN_PER_TET + nbConstraints + coord) = 1;

    knowns(nbTets * NB_UNKNOWN_PER_TET + nbConstraints) = 0;
    knowns(nbTets * NB_UNKNOWN_PER_TET + nbConstraints + 1) = 0;
  }

  // std::cout << system << std::endl;
  // std::cout << knowns << std::endl;
  system.makeCompressed();

  Eigen::ConjugateGradient<Eigen::SparseMatrix<double>,
                           Eigen::Lower | Eigen::Upper>
      solver;

  solver.setMaxIterations(1000);
  solver.setTolerance(1e-8);
  solver.compute(system);

  Eigen::VectorXd result;
  if (solver.info() != Eigen::Success) {
    std::cerr << "Failed to initialize CG solver\n";
  } else {
    result = solver.solve(knowns);

    std::cerr << "Iterations: " << solver.iterations() << '\n';
    std::cerr << "Estimated error: " << solver.error() << '\n';

    if (solver.info() != Eigen::Success) {
      std::cerr << "CG failed to converge\n";
    }
  }

  size_t dimensionOffset = NB_UNKNOWN_PER_TET;
  auto setLocalVertices = [&](size_t tetId, size_t localNodeId, Vertex value) {
    size_t offset = tetId * NB_UNKNOWN_PER_TET;
    for (size_t dim = 0; dim < DIMENSION; dim++)
      result[offset + dim * NB_NODE_PER_TET + localNodeId] = value[dim];
  };
  auto getLocalVertices = [&](size_t tetId, size_t localNodeId) -> Vertex {
    size_t offset = tetId * NB_UNKNOWN_PER_TET;
    Vertex res;
    for (size_t dim = 0; dim < DIMENSION; dim++)
      res[dim] = result[offset + dim * NB_NODE_PER_TET + localNodeId];
    return res;
  };
  auto getGlobalVertices = [&](size_t tetId, size_t globalNodeId) {
    auto tetNodes = mesh.tets[tetId];
    for (size_t i = 0; i < tetNodes.size(); i++) {
      if (tetNodes[i] == globalNodeId)
        return getLocalVertices(tetId, i);
    }
    assert(false);
  };
  auto setGlobalVertices = [&](size_t tetId, size_t globalNodeId,
                               Vertex value) {
    auto tetNodes = mesh.tets[tetId];
    for (size_t i = 0; i < tetNodes.size(); i++) {
      if (tetNodes[i] == globalNodeId)
        return setLocalVertices(tetId, i, value);
    }
    assert(false);
  };

  std::vector<Eigen::Vector3d> vertices(mesh.tets.size() * NB_NODE_PER_TET);
  for (size_t i = 0; i < mesh.tets.size(); i++) {
    for (size_t j = 0; j < NB_NODE_PER_TET; j++)
      vertices[NB_NODE_PER_TET * i + j] = getLocalVertices(i, j);
  }
  saveMeshVtk(vertices, "nospanning.vtk");
  // std::cout << "RESULT " << std::endl << result << std::endl;

  std::vector<bool> include(mesh.tets.size());
  std::fill(include.begin(), include.end(), false);
  for (size_t faceId : spanTree) {
    Face &face = mesh.faces[faceId];
    auto faceNodes = face.vertices;
    if (face.adjacentCells[1] == UNDEFINED_ELEMENT)
      continue;
    size_t parent = face.adjacentCells[0];
    size_t child = face.adjacentCells[1];
    include[parent] = true;
    include[child] = true;

    // Find the child tip
    auto childNodes = mesh.tets[child];
    size_t childTipLocalId = NB_NODE_PER_TET;
    for (childTipLocalId = 0; childTipLocalId < NB_NODE_PER_TET;
         childTipLocalId++) {
      size_t tipGlobalId = childNodes[childTipLocalId];
      if (std::find(faceNodes.begin(), faceNodes.end(), tipGlobalId) ==
          faceNodes.end())
        break;
    }
    assert(childTipLocalId < NB_NODE_PER_TET);
    Vertex childTip(getLocalVertices(child, childTipLocalId));

    Vertex oldFace[NB_NODE_PER_FACE];
    Vertex newFace[NB_NODE_PER_FACE];

    for (size_t i = 0; i < NB_NODE_PER_FACE; i++) {
      size_t node = faceNodes[i];
      Vertex fp = getGlobalVertices(parent, node);
      Vertex fc = getGlobalVertices(child, node);
      oldFace[i] = fc;
      newFace[i] = fp;
      setGlobalVertices(child, node, fp);
    }

    auto newTip = tipProject(oldFace, newFace, childTip);
    setLocalVertices(child, childTipLocalId, newTip);
  }

  vertices.clear();
  for (size_t i = 0; i < mesh.tets.size(); i++) {
    if (!include[i])
      continue;
    for (size_t j = 0; j < NB_NODE_PER_TET; j++)
      vertices.push_back(getLocalVertices(i, j));
  }
  saveMeshGmsh(vertices);
  saveMeshVtk(vertices);

  return EXIT_SUCCESS;
}
