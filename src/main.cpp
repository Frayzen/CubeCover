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

static double computeArea(TetrahedralMesh &mesh, size_t tetId) {
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
  double data[3][3] = {-999};
  Tetrahedron &currentTet = mesh.tets[tetId];
  Vertex &p0 = mesh.vertices[currentTet[0]];
  for (size_t i = 1; i < NB_NODE_PER_TET; i++) {
    Vertex &pi = mesh.vertices[currentTet[i]];
    for (size_t coord = 0; coord < DIMENSION; coord++)
      data[i - 1][coord] = pi[coord] - p0[coord];
  }
  Eigen::Map<Eigen::Matrix<double, 3, 3, Eigen::ColMajor>> mat(data[0]);
  auto invT = mat.inverse();

  DiffTetShapeFunctions res{
      .g0 = {},
      .g1 = invT.row(0),
      .g2 = invT.row(1),
      .g3 = invT.row(2),
  };
  res.g0 = -(res.g1 + res.g2 + res.g3);

  // std::cout << "=============================" << std::endl;
  // std::cout << "GOT mat " << std::endl << mat << std::endl;
  // std::cout << "GOT AREA " << computeArea(mesh, tetId) << std::endl;
  // std::cout << "g0 " << std::endl << res.g0 << std::endl << std::endl;
  // std::cout << "g1 " << std::endl << res.g1 << std::endl << std::endl;
  // std::cout << "g2 " << std::endl << res.g2 << std::endl << std::endl;
  // std::cout << "g3 " << std::endl << res.g3 << std::endl << std::endl;
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
  convertToSpanTree(mesh);

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
  size_t nbConstraints =
      mesh.amountInternalFace * (NB_NODE_PER_FACE + 1) * DIMENSION;

  size_t bcAmount = DIMENSION;
  size_t N = nbTets * NB_UNKNOWN_PER_TET + nbConstraints + bcAmount;

  // The unknowns are in the shape [u1 u2 ... v1 v2 ... w1 w2 ... u5 u6 ... v5
  // v6 ... w5 w6 ... deltas ... ]
  Eigen::SparseMatrix<double> system(N, N);
  system.reserve(nbTets * NB_UNKNOWN_PER_TET + nbConstraints);

  Eigen::VectorX<double> knowns(N);

  // Laplace matrix

  for (size_t i = 0; i < mesh.tets.size(); i++) {
    DiffTetShapeFunctions shapeFuncs = computeDiffShapeFunctions(mesh, i);
    Eigen::Matrix3d frame;   // TODO set the frame not constant
    frame.col(0) << 1, 0, 0; // U
    frame.col(1) << 0, 1, 0; // V
    frame.col(2) << 0, 0, 1; // W

    for (size_t nodeId = 0; nodeId < NB_NODE_PER_FACE; nodeId++) {
      size_t base = i * NB_UNKNOWN_PER_TET + nodeId * NB_NODE_PER_TET;
      for (size_t di = 0; di < NB_NODE_PER_TET; di++) {
        for (size_t dj = 0; dj < NB_NODE_PER_TET; dj++) {
          double coeff = shapeFuncs[di].transpose() * shapeFuncs[dj];
          system.insert(base + di, base + dj) = coeff;
          knowns(base + di) = shapeFuncs[di].transpose() * frame.col(nodeId);
        }
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

    // Here, we assume that Pi is Identity
    for (size_t coord = 0; coord < DIMENSION; coord++) {
      size_t colTetT =
          face.adjacentCells[0] * NB_UNKNOWN_PER_TET + coord * NB_NODE_PER_TET;
      size_t colTetS =
          face.adjacentCells[1] * NB_UNKNOWN_PER_TET + coord * NB_NODE_PER_TET;
      for (size_t nodeId = 0; nodeId < NB_NODE_PER_FACE; nodeId++) {
        size_t row = constraintOffset + constraintId;
        constraintId++;

        // Enforce f(u_t) - f(v_t) - Pi(f(u_s) - f(v_s)) = 0
        knowns(row) = 0;
        size_t nextId = (nodeId + 1) % 4;
        // C
        system.insert(row, colTetT + nodeId) = 1;
        system.insert(row, colTetS + nodeId) = -1; // TODO apply Pi
        system.insert(row, colTetT + nextId) = -1;
        system.insert(row, colTetS + nextId) = 1; // TODO apply Pi
        // C^t
        system.insert(colTetT + nodeId, row) = 1;
        system.insert(colTetS + nodeId, row) = -1; // TODO apply Pi
        system.insert(colTetT + nextId, row) = -1;
        system.insert(colTetS + nextId, row) = 1; // TODO apply Pi
      }

      // Enforce f(u_t)
      size_t row = constraintOffset + constraintId;
      constraintId++;
      knowns(row) = 0;
      // C
      system.insert(row, colTetT) = 1;
      system.insert(row, colTetS) = -1; // TODO apply Pi
      // C^t
      system.insert(colTetT, row) = 1;
      system.insert(colTetS, row) = -1; // TODO apply Pi
    }
  }
  assert(constraintId == nbConstraints);

  // Boundary conditions
  for (size_t coord = 0; coord < DIMENSION; coord++) {
    system.insert(nbTets * NB_UNKNOWN_PER_TET + nbConstraints + coord,
                  coord * NB_NODE_PER_TET) = 1;
    // C^t
    system.insert(coord * NB_NODE_PER_TET,
                  nbTets * NB_UNKNOWN_PER_TET + nbConstraints + coord) = 1;
    knowns(nbTets * NB_UNKNOWN_PER_TET + nbConstraints) = 0;
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

  for (size_t faceId = 0; faceId < mesh.faces.size(); faceId++) {
    Face &face = mesh.faces[faceId];
    auto faceNodes = face.vertices;
    if (face.adjacentCells[1] == UNDEFINED_ELEMENT)
      continue;
    size_t parent = face.adjacentCells[0];
    size_t child = face.adjacentCells[1];

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

    double distToTip[NB_NODE_PER_FACE];
    Vertex faceNodePos[NB_NODE_PER_FACE];

    for (size_t i = 0; i < NB_NODE_PER_FACE; i++) {
      size_t node = faceNodes[i];
      Vertex fp = getGlobalVertices(parent, node);
      Vertex fc = getGlobalVertices(child, node);
      distToTip[i] = (fc - childTip).norm();
      setGlobalVertices(child, node, fp);
      faceNodePos[i] = fp;
    }

    auto [newTipPos1, newTipPos2] = trilaterate3D(faceNodePos, distToTip);
    setLocalVertices(child, childTipLocalId, newTipPos1);
  }

  std::vector<Eigen::Vector3d> vertices(mesh.tets.size() * NB_NODE_PER_TET);
  for (size_t i = 0; i < mesh.tets.size(); i++) {
    for (size_t j = 0; j < NB_NODE_PER_TET; j++)
      vertices[NB_NODE_PER_TET * i + j] = getLocalVertices(i, j);
  }
  saveMesh(vertices);

  return EXIT_SUCCESS;
}
