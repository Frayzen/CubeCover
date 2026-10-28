
#include "iterator.hpp"
#include "constants.hpp"
#include "mesh_saver.hpp"
#include "shape_functions.hpp"
#include "vertex.hpp"
#include <Eigen/IterativeLinearSolvers>

size_t CubeCoverIterator::getId(size_t tetId, size_t localNode, size_t coord) {
  const size_t vertexId = mesh.tets[tetId][localNode];
  return vertexId * DIMENSION + coord;
}

void CubeCoverIterator::fillLaplaceTop() {
  for (size_t i = 0; i < mesh.tets.size(); ++i) {
    const DiffTetShapeFunctions shapeFuncs = computeDiffShapeFunctions(mesh, i);

    const double vol = computeVolume(mesh, i);
    const Eigen::Matrix3d frame = frames[i];

    for (size_t coord = 0; coord < DIMENSION; ++coord) {
      for (size_t di = 0; di < NB_NODE_PER_TET; ++di) {
        const size_t row = getId(i, di, coord);

        for (size_t dj = 0; dj < NB_NODE_PER_TET; ++dj) {
          const size_t col = getId(i, dj, coord);

          const double coeff = vol * shapeFuncs[di].dot(shapeFuncs[dj]);

          A.coeffRef(row, col) += coeff;
        }

        b(row) += vol * shapeFuncs[di].dot(frame.col(coord));
      }
    }
  }
}

void CubeCoverIterator::fillCornerConstraints() {
  size_t offset = sizeLaplace;
  size_t constraintId = 0;

  // std::array<std::array<size_t, 4>, 6> sharingCoord = {};
  // sharingCoord[0] = {0, 1, 2, 3}; // X+
  // sharingCoord[1] = {4, 5, 6, 7}; // X-
  // sharingCoord[2] = {0, 1, 4, 5}; // Y+
  // sharingCoord[3] = {2, 3, 6, 7}; // Y-
  // sharingCoord[4] = {0, 3, 4, 7}; // Z+
  // sharingCoord[5] = {1, 2, 5, 6}; // Z-

  // for (size_t i = 0; i < 2; i++) {
  //   size_t d = i / 2;
  //   size_t colFirstNode = cornerIds[sharingCoord[i][0]] * DIMENSION + d;
  //   for (size_t j = 1; j < 4; j++) {
  //     size_t col = cornerIds[sharingCoord[i][j]] * DIMENSION + d;
  //     size_t row = offset + constraintId;
  //     A.coeffRef(row, col) = 1;
  //     A.coeffRef(col, row) = 1;

  //     A.coeffRef(row, colFirstNode) = -1;
  //     A.coeffRef(colFirstNode, row) = -1;

  //     b.coeffRef(row) = 0;
  //     constraintId++;
  //   }
  // }

  for (size_t d = 0; d < DIMENSION; d++) {
    size_t row = offset + constraintId;
    size_t col = cornerIds[0] * DIMENSION + d;
    A.coeffRef(row, col) = 1;
    A.coeffRef(col, row) = 1;
    b.coeffRef(row) = 0;
    constraintId++;
  }

  assert(constraintId == sizeConstraint);
}

void CubeCoverIterator::initialize() {
  std::cout << "Initializing...\n";
  fillLaplaceTop();
  fillCornerConstraints();
  A.makeCompressed();
}
void CubeCoverIterator::iterate() {
  solve();
  saveMesh("nospanning.vtu");
}

void CubeCoverIterator::solve() {
  std::cout << "Solving..." << std::endl;

  Eigen::LeastSquaresConjugateGradient<Eigen::SparseMatrix<double>> solver;

  solver.setMaxIterations(10000);
  solver.setTolerance(1e-15);

  solver.compute(A);

  if (solver.info() != Eigen::Success) {
    std::cerr << "Failed to initialize solver\n";
    throw std::runtime_error("Solver initialization failed");
  }

  x = solver.solveWithGuess(b, x);

  std::cerr << "Iterations: " << solver.iterations() << '\n';
  std::cerr << "Estimated error: " << solver.error() << '\n';

  const double residual = (A * x - b).norm();

  std::cerr << "Residual: " << residual << '\n';

  if (solver.info() != Eigen::Success) {
    std::cerr << "Least squares solver failed\n";
    throw std::runtime_error("Least squares solver failed");
  }
}

Vertex CubeCoverIterator::getVertex(size_t tetId, size_t localNodeId) {
  Vertex res;
  for (size_t i = 0; i < DIMENSION; i++) {
    res[i] = x[getId(tetId, localNodeId, i)];
  }
  return res;
}

std::array<Vertex, NB_NODE_PER_TET>
CubeCoverIterator::getTetrahedraVertices(size_t tetId) {
  std::array<Vertex, NB_NODE_PER_TET> res;
  for (size_t i = 0; i < NB_NODE_PER_TET; i++)
    res[i] = getVertex(tetId, i);
  return res;
}

void CubeCoverIterator::saveMesh(std::string path) {
  std::vector<Vertex> vertices;
  Vertex val;
  for (size_t i = 0; i < mesh.vertices.size(); i++) {
    size_t base = DIMENSION * i;
    val << x[base], x[base + 1], x[base + 2];
    vertices.push_back(val);
  }
  saveMeshVtk(vertices, mesh.tets, frames, path);
}
