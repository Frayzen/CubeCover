
#include "iterator.hpp"
#include "frame_field.hpp"
#include "shape_functions.hpp"
#include <Eigen/IterativeLinearSolvers>

size_t CubeCoverIterator::getRow(size_t tetId, size_t localNode, size_t coord) {
  return tetId * NB_UNKNOWN_PER_TET + localNode * DIMENSION + coord;
}

void CubeCoverIterator::fillLaplaceTop() {
  for (size_t i = 0; i < mesh.tets.size(); i++) {
    DiffTetShapeFunctions shapeFuncs = computeDiffShapeFunctions(mesh, i);
    Eigen::Matrix3d frame = frames[i];

    const double vol = computeVolume(mesh, i);

    for (size_t coord = 0; coord < DIMENSION; ++coord) {
      for (size_t di = 0; di < NB_NODE_PER_TET; ++di) {
        const size_t row = getRow(i, di, coord);

        for (size_t dj = 0; dj < NB_NODE_PER_TET; ++dj) {
          const size_t col = getRow(i, dj, coord);

          const double coeff = vol * shapeFuncs[di].dot(shapeFuncs[dj]);

          A.coeffRef(row, col) += coeff;
        }

        b(row) += vol * shapeFuncs[di].dot(frame.col(coord));
      }
    }
  }
}

void CubeCoverIterator::fillConstraintTop() {
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

  const size_t constraintOffset = mesh.tets.size() * NB_UNKNOWN_PER_TET;
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
      auto [coordS, sign] = findCoord(pi, coord);

      for (size_t nodeId = 0; nodeId < NB_NODE_PER_FACE; nodeId++) {
        size_t node = face.vertices[nodeId];
        size_t nextId = (nodeId + 1) % NB_NODE_PER_FACE;
        size_t next = face.vertices[nextId];

        size_t row = constraintOffset + constraintId;
        constraintId++;

        // Enforce f(u_t) - f(v_t) - Pi(f(u_s) - f(v_s)) = 0
        b(row) = 0;

        // Node-major indexing for tetrahedron T
        size_t colTetTNode = getRow(tetTId, tetT.getLocalId(node), coord);
        size_t colTetTNext = getRow(tetTId, tetT.getLocalId(next), coord);

        // Node-major indexing for tetrahedron S
        size_t colTetSNode = getRow(tetSId, tetS.getLocalId(node), coordS);
        size_t colTetSNext = getRow(tetSId, tetS.getLocalId(next), coordS);

        // C
        A.insert(row, colTetTNode) = 1;
        A.insert(row, colTetSNode) = -sign;
        A.insert(row, colTetTNext) = -1;
        A.insert(row, colTetSNext) = sign;

        // C^t
        A.insert(colTetTNode, row) = 1;
        A.insert(colTetSNode, row) = -sign;
        A.insert(colTetTNext, row) = -1;
        A.insert(colTetSNext, row) = sign;
      }
    }
  }

  assert(constraintId == sizeConstraint);
}

void CubeCoverIterator::initialize() {
  std::cout << "Initializing...\n";
  fillLaplaceTop();
  fillConstraintTop();
  A.makeCompressed();
}

void CubeCoverIterator::solve() {
  std::cout << "Solving..." << std::endl;

  Eigen::ConjugateGradient<Eigen::SparseMatrix<double>,
                           Eigen::Lower | Eigen::Upper>
      solver;

  solver.setMaxIterations(10000);
  solver.setTolerance(1e-4);
  solver.compute(A);

  Eigen::VectorXd result;

  if (solver.info() != Eigen::Success) {
    std::cerr << "Failed to initialize CG solver\n";
  } else {
    x = solver.solveWithGuess(b, x);

    std::cerr << "Iterations: " << solver.iterations() << '\n';
    std::cerr << "Estimated error: " << solver.error() << '\n';

    if (solver.info() != Eigen::Success) {
      std::cerr << "CG failed to converge\n";
    }
  }
}
