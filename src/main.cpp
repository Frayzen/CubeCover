#include "constants.hpp"
#include "mesh_loader.hpp"
#include <Eigen/Core>
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
  auto inv = mat.inverse();

  DiffTetShapeFunctions res{
      .g0 = {},
      .g1 = inv.row(0),
      .g2 = inv.row(1),
      .g3 = inv.row(2),
  };
  res.g0 = Vertex(1, 1, 1) - res.g1 - res.g2 - res.g3;

  // std::cout << "=============================" << std::endl;
  // std::cout << "GOT AREA " << computeArea(mesh, tetId) << std::endl;
  // std::cout << res.g0 << std::endl;
  // std::cout << res.g1 << std::endl;
  // std::cout << res.g2 << std::endl;
  // std::cout << res.g3 << std::endl;
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
  size_t nbConstraints = mesh.amountInternalFace * NB_NODE_PER_FACE;

  // the one more is for the constraint that sets one vertex coordinate
  size_t N = nbTets * NB_UNKNOWN_PER_TET + nbConstraints + 1;
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
          system.coeffRef(base + di, base + dj) = coeff;
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
    if (face.adjacentCells[1] == UNDEFINED_ADJACENT_CELL)
      continue;
    // Here, we assume that Pi is Identity
    for (size_t nodeId = 0; nodeId < NB_NODE_PER_FACE; nodeId++) {
      size_t row = constraintOffset + constraintId;
      constraintId++;
      size_t colTetT = face.adjacentCells[0] * NB_UNKNOWN_PER_TET + nodeId;
      size_t colTetS = face.adjacentCells[1] * NB_UNKNOWN_PER_TET + nodeId;

      knowns(row) = 0;

      // C
      system.coeffRef(row, colTetT) = 1;
      system.coeffRef(row, colTetS) = -1; // TODO apply Pi
      // C^t
      system.coeffRef(colTetT + nodeId, row) = 1;
      system.coeffRef(colTetS + nodeId, row) = -1; // TODO apply Pi

      if (nodeId == 0) { // Enforce f(u_t) - Pi(f(u_s)) = 0
        continue;
      } else { // Enforce f(u_t) - f(v_t) - Pi(f(u_s) - f(v_s)) = 0
        system.coeffRef(row, colTetT - 1) = -1;
        system.coeffRef(row, colTetS - 1) = 1; // TODO apply Pi
        // C^t
        system.coeffRef(colTetT - 1, row) = -1;
        system.coeffRef(colTetS - 1, row) = 1; // TODO apply Pi
      }
    }
  }
  assert(constraintId == nbConstraints);

  // Boundary conditions
  system.coeffRef(nbTets * NB_UNKNOWN_PER_TET + nbConstraints, 0) = 1;
  knowns(nbTets * NB_UNKNOWN_PER_TET + nbConstraints) = 0;

  system.makeCompressed();

  // Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
  // solver.analyzePattern(system);
  // solver.factorize(system);
  // if (solver.info() != Eigen::Success) {
  //   std::cout << "Cannot solve" << "\n";
  // }
  // Eigen::VectorXd solution = solver.solve(knowns);

  return EXIT_SUCCESS;
}
