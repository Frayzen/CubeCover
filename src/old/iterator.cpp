#include "iterator.hpp"
#include "constants.hpp"
#include "mesh_saver.hpp"
#include "shape_functions.hpp"
#include "signed_axis.hpp"
#include "ids.hpp"
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SparseLU>

size_t CubeCoverIterator::getId(size_t tetId, size_t localNode, size_t coord) {
  const size_t vertexId = mesh.tets[tetId][localNode];
  return vertexId * DIMENSION + coord;
}

void CubeCoverIterator::fillLaplace() {
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

void CubeCoverIterator::fillConstraints() {
  size_t offset = sizeLaplace;
  size_t constraintId = 0;

  for (NodeId nodeId = 0; nodeId < mesh.nbBoundaryNode; nodeId++) {
    for (Axis d = 0; d < DIMENSION; d++) {
      if (nodeToAxisConstraints[nodeId][d] == ZERO)
        continue;
      size_t row = offset + constraintId;
      size_t col = nodeId * DIMENSION + d;
      A.coeffRef(row, col) = 1;
      A.coeffRef(col, row) = 1;
      b.coeffRef(row) = getIntegerSignValue(nodeToAxisConstraints[nodeId][d]);
      constraintId++;
    }
  }

  assert(constraintId == sizeConstraint);
}

void CubeCoverIterator::initialize() {
  std::cout << "Initializing...\n";
  fillLaplace();
  fillConstraints();
  A.makeCompressed();
}
void CubeCoverIterator::solve() {
  std::cout << "Solving..." << std::endl;

  // Eigen::LeastSquaresConjugateGradient<Eigen::SparseMatrix<double>> solver;

  // solver.setMaxIterations(10000);
  // solver.setTolerance(1e-15);

  // solver.compute(A);

  // if (solver.info() != Eigen::Success) {
  //   std::cerr << "Failed to initialize solver\n";
  //   throw std::runtime_error("Solver initialization failed");
  // }

  // x = solver.solveWithGuess(b, x);

  // std::cerr << "Iterations: " << solver.iterations() << '\n';
  // std::cerr << "Estimated error: " << solver.error() << '\n';

  // const double residual = (A * x - b).norm();

  // std::cerr << "Residual: " << residual << '\n';

  // if (solver.info() != Eigen::Success) {
  //   std::cerr << "Least squares solver failed\n";
  //   throw std::runtime_error("Least squares solver failed");
  // }

  Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;

  solver.compute(A);

  if (solver.info() != Eigen::Success) {
    throw std::runtime_error("Factorization failed");
  }

  x = solver.solve(b);

  if (solver.info() != Eigen::Success) {
    throw std::runtime_error("Solving failed");
  }

  double maxViolation = 0.0;

  for (NodeId nodeId = 0; nodeId < mesh.nbBoundaryNode; nodeId++) {
    for (Axis d = 0; d < DIMENSION; d++) {
      auto constraint = nodeToAxisConstraints[nodeId][d];

      if (constraint == ZERO)
        continue;

      size_t col = nodeId * DIMENSION + d;

      double expected = getIntegerSignValue(constraint);

      maxViolation = std::max(maxViolation, std::abs(x[col] - expected));
    }
  }

  std::cout << "Maximum constraint violation: " << maxViolation << '\n';
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
  std::vector<Eigen::Vector3d> frameX(mesh.tets.size());
  std::vector<Eigen::Vector3d> frameY(mesh.tets.size());
  std::vector<Eigen::Vector3d> frameZ(mesh.tets.size());

  std::array<std::vector<int>, DIMENSION> constraints;
  for (Axis d = 0; d < DIMENSION; d++)
    constraints[d] = std::vector<int>(mesh.vertices.size(), 0);

  for (NodeId nodeId = 0; nodeId < mesh.nbBoundaryNode; nodeId++) {
    for (Axis d = 0; d < DIMENSION; d++) {
      if (nodeToAxisConstraints[nodeId][d] == ZERO)
        continue;
      constraints[d][nodeId] =
          getIntegerSignValue(nodeToAxisConstraints[nodeId][d]);
    }
  }

  for (size_t tetId = 0; tetId < mesh.tets.size(); tetId++) {
    frameX[tetId] = frames[tetId].col(0);
    frameY[tetId] = frames[tetId].col(1);
    frameZ[tetId] = frames[tetId].col(2);
  }

  TetrahedralMesh newMesh(mesh);
  newMesh.vertices = vertices;
  MeshSaver saver(newMesh, path);
  saver.addField<ElementType::CellType>("frameX", frameX);
  saver.addField<ElementType::CellType>("frameY", frameY);
  saver.addField<ElementType::CellType>("frameZ", frameZ);

  saver.addField<ElementType::NodeType>("constraintX", constraints[0]);
  saver.addField<ElementType::NodeType>("constraintY", constraints[1]);
  saver.addField<ElementType::NodeType>("constraintZ", constraints[2]);

  saver.saveMesh();
}
