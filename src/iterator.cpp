
#include "iterator.hpp"
#include "frame_field.hpp"
#include "mesh_saver.hpp"
#include "shape_functions.hpp"
#include "span_tree.hpp"
#include "trilaterate.hpp"
#include <Eigen/IterativeLinearSolvers>
#include <system_error>

// Deduplicate DOFs by global mesh vertex.
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

void CubeCoverIterator::fillConstraintTop() {
  auto findCoord = [](Eigen::Matrix3i &pi,
                      size_t coord) -> std::tuple<size_t, int> {
    Eigen::Vector3i direction = Eigen::Vector3i::Zero();
    direction[coord] = 1;
    direction = pi * direction;

    for (size_t i = 0; i < DIMENSION; ++i) {
      if (direction[i] != 0)
        return {i, direction[i]};
    }

    std::cerr << "Error in PI " << pi << std::endl;
    assert(false);
    return {0, 0};
  };

  const size_t constraintOffset = mesh.vertices.size() * DIMENSION;

  size_t constraintId = 0;

  for (size_t faceId = 0; faceId < mesh.faces.size(); ++faceId) {
    Face &face = mesh.faces[faceId];

    if (face.adjacentCells[1] == UNDEFINED_ELEMENT)
      continue;

    const size_t tetTId = face.adjacentCells[0];
    const Tetrahedron &tetT = mesh.tets[tetTId];

    const size_t tetSId = face.adjacentCells[1];
    const Tetrahedron &tetS = mesh.tets[tetSId];

    const Eigen::Matrix3i pi = computePiST(frames[tetTId], frames[tetSId]);

    for (size_t coord = 0; coord < DIMENSION; ++coord) {
      const auto [coordS, sign] =
          findCoord(const_cast<Eigen::Matrix3i &>(pi), coord);

      for (size_t nodeId = 0; nodeId < NB_NODE_PER_FACE; ++nodeId) {
        const size_t node = face.vertices[nodeId];
        const size_t nextId = (nodeId + 1) % NB_NODE_PER_FACE;
        const size_t next = face.vertices[nextId];

        const size_t row = constraintOffset + constraintId++;

        b(row) = 0.0;

        const size_t colTetTNode = getId(tetTId, tetT.getLocalId(node), coord);
        const size_t colTetTNext = getId(tetTId, tetT.getLocalId(next), coord);

        const size_t colTetSNode = getId(tetSId, tetS.getLocalId(node), coordS);
        const size_t colTetSNext = getId(tetSId, tetS.getLocalId(next), coordS);

        // Constraint matrix C.
        A.coeffRef(row, colTetTNode) += 1.0;
        A.coeffRef(row, colTetSNode) -= sign;
        A.coeffRef(row, colTetTNext) -= 1.0;
        A.coeffRef(row, colTetSNext) += sign;

        // Transpose C^T.
        A.coeffRef(colTetTNode, row) += 1.0;
        A.coeffRef(colTetSNode, row) -= sign;
        A.coeffRef(colTetTNext, row) -= 1.0;
        A.coeffRef(colTetSNext, row) += sign;
      }
    }
  }

  assert(constraintId == sizeConstraint);
}

void CubeCoverIterator::buildCoordsToRemove() {
  coordsToRemove.clear();
  std::vector<size_t> toTreat;

  for (const Face &face : mesh.faces) {
    if (face.adjacentCells[1] != UNDEFINED_ELEMENT)
      continue;
    toTreat.push_back(face.adjacentCells[0]);
  }

  std::sort(toTreat.begin(), toTreat.end(),
            [&](const size_t &aId, const size_t &bId) {
              const Tetrahedron &a = mesh.tets[aId];
              const Tetrahedron &b = mesh.tets[bId];
              const auto aBary = a.getBarycenter(mesh.vertices);
              const auto bBary = b.getBarycenter(mesh.vertices);
              return aBary.norm() < bBary.norm();
            });

  for (size_t tetId : toTreat) {
    const Tetrahedron &tet = mesh.tets[tetId];
    const auto &facesId = mesh.tetToFaces[tetId];

    // One boolean per node and coordinate dimension.
    bool toAdd[NB_NODE_PER_TET][DIMENSION] = {};

    for (size_t faceId : facesId) {
      const Face &face = mesh.faces[faceId];

      Vertex nodes[NB_NODE_PER_FACE];
      for (size_t i = 0; i < NB_NODE_PER_FACE; ++i)
        nodes[i] = mesh.vertices[face.vertices[i]];

      Vertex normal = getNormal(nodes);
      normal.normalize();

      // Find the two coordinate axes least aligned with the normal.
      std::array<size_t, DIMENSION> dims = {0, 1, 2};
      std::sort(dims.begin(), dims.end(), [&](size_t a, size_t b) {
        return std::abs(normal[a]) < std::abs(normal[b]);
      });

      for (size_t i = 0; i < NB_NODE_PER_FACE; ++i) {
        const size_t localNode = tet.getLocalId(face.vertices[i]);

        for (size_t j = 0; j < 2; ++j)
          toAdd[localNode][dims[j]] = true;
      }
    }

    for (size_t i = 0; i < NB_NODE_PER_TET; ++i) {
      for (size_t dim = 0; dim < DIMENSION; ++dim) {
        if (toAdd[i][dim]) {
          coordsToRemove.push_back(getId(tetId, i, dim));
        }
      }
    }
  }
}

void CubeCoverIterator::initialize() {
  std::cout << "Initializing...\n";
  fillLaplaceTop();
  fillConstraintTop();
  A.makeCompressed();
  spanTree = convertToSpanTree(mesh);
  buildCoordsToRemove();
}
void CubeCoverIterator::iterate() {
  solve();
  saveMesh("nospanning.vtu");
  applySpanTree();
  removeOneCoord(coordsToRemove[coordsToRemove.size() - 1]);
  coordsToRemove.pop_back();
}

void CubeCoverIterator::solve() {
  std::cout << "Solving..." << std::endl;

  Eigen::ConjugateGradient<Eigen::SparseMatrix<double>,
                           Eigen::Lower | Eigen::Upper>
      solver;

  solver.setMaxIterations(10000);
  solver.setTolerance(1e-5);
  solver.compute(A);

  Eigen::VectorXd result;

  if (solver.info() != Eigen::Success) {
    std::cerr << "Failed to initialize CG solver\n";
  } else {
    x = solver.solveWithGuess(b, x);
    Eigen::VectorXd r = A * x - b;

    std::cerr << "Iterations: " << solver.iterations() << '\n';
    std::cerr << "Estimated error: " << solver.error() << '\n';

    if (solver.info() != Eigen::Success) {
      std::cerr << "CG failed to converge\n";
      throw std::system_error();
    }
  }
}

Vertex CubeCoverIterator::getVertex(size_t tetId, size_t localNodeId) {
  Vertex res;
  for (size_t i = 0; i < DIMENSION; i++) {
    size_t id = getId(tetId, localNodeId, i);
    bool inFixedValues = scalarLocation[id] < 0;
    if (inFixedValues)
      res[i] = fixedValues[-scalarLocation[id] - 1];
    else
      res[i] = x[scalarLocation[id]];
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

void CubeCoverIterator::applySpanTree() {
  Eigen::VectorX<double> result = x;
  auto setVertex = [&](size_t tetId, size_t localNodeId, const Vertex &value) {
    for (size_t dim = 0; dim < DIMENSION; ++dim) {
      const size_t id = getId(tetId, localNodeId, dim);
      const auto loc = scalarLocation[id];

      // Fixed coordinates must never be modified.
      if (loc >= 0)
        result[loc] = value[dim];
    }
  };

  auto getVertexFromResult = [&](size_t tetId, size_t localNodeId) -> Vertex {
    Vertex res;

    for (size_t dim = 0; dim < DIMENSION; ++dim) {
      const size_t id = getId(tetId, localNodeId, dim);
      const auto loc = scalarLocation[id];

      res[dim] = (loc < 0) ? fixedValues[-loc - 1] : result[loc];
    }

    return res;
  };

  for (size_t faceId : spanTree) {
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
    size_t parentTipLocalId = NB_NODE_PER_TET;
    auto parentNodes = mesh.tets[parent];

    for (parentTipLocalId = 0; parentTipLocalId < NB_NODE_PER_TET;
         parentTipLocalId++) {
      size_t tipGlobalId = parentNodes[parentTipLocalId];

      if (std::find(faceNodes.begin(), faceNodes.end(), tipGlobalId) ==
          faceNodes.end())
        break;
    }

    assert(parentTipLocalId < NB_NODE_PER_TET);

    Vertex childTip(getVertex(child, childTipLocalId));

    Vertex oldFace[NB_NODE_PER_FACE];
    Vertex newFace[NB_NODE_PER_FACE];

    for (size_t i = 0; i < NB_NODE_PER_FACE; i++) {
      size_t node = faceNodes[i];

      Vertex fp = getVertexFromResult(parent, parentNodes.getLocalId(node));
      Vertex fc = getVertex(child, childNodes.getLocalId(node));

      oldFace[i] = fc;
      newFace[i] = fp;

      for (size_t i = 0; i < DIMENSION; i++)
        setVertex(child, childNodes.getLocalId(node), fp);
    }

    Vertex &parentTip = mesh.vertices[parentNodes[parentTipLocalId]];
    try {
      auto newTip = tipProject(oldFace, newFace, childTip);
      for (size_t i = 0; i < DIMENSION; i++)
        setVertex(child, childTipLocalId, newTip);
    } catch (std::exception &e) {
      std::cerr << "Error with face " << faceId << std::endl;
      std::cerr << "Parent tet " << parent << std::endl;
      std::cerr << "Child tet " << child << std::endl;
      throw e;
    }
  }
  x = result;
}

void CubeCoverIterator::saveMesh(std::string path) {
  std::vector<Vertex> vertices;
  std::vector<Eigen::Vector3i> fixedMasks;

  vertices.reserve(mesh.tets.size() * NB_NODE_PER_TET);
  fixedMasks.reserve(mesh.tets.size() * NB_NODE_PER_TET);

  for (size_t tetId = 0; tetId < mesh.tets.size(); ++tetId) {
    for (size_t localNode = 0; localNode < NB_NODE_PER_TET; ++localNode) {
      vertices.push_back(getVertex(tetId, localNode));

      Eigen::Vector3i mask;
      for (size_t dim = 0; dim < DIMENSION; ++dim) {
        const size_t id = getId(tetId, localNode, dim);
        mask[dim] = (scalarLocation[id] < 0) ? 1 : 0;
      }

      fixedMasks.push_back(mask);
    }
  }

  saveMeshVtk(vertices, frames, fixedMasks, path);
}

// void CubeCoverIterator::saveMesh(std::string path) {
//   std::vector<Vertex> vertices;
//   vertices.reserve(mesh.tets.size() * NB_NODE_PER_TET);
//   for (size_t i = 0; i < mesh.tets.size(); i++) {
//     auto toAdd = getTetrahedraVertices(i);
//     vertices.insert(vertices.end(), toAdd.begin(), toAdd.end());
//   }
//   saveMeshVtk(vertices, frames, path);
// }

void CubeCoverIterator::removeOneCoord(size_t coordId) {
  const Eigen::Index k = scalarLocation[coordId];

  // This coordinate has already been fixed.
  if (k < 0)
    return;

  const double value = std::round(x[k]);
  const Eigen::Index n = A.rows();

  assert(A.cols() == n);
  assert(b.size() == n);
  assert(k < n);

  // Eliminate the fixed variable's contribution from b.
  Eigen::VectorXd newB = b - A.col(k) * value;

  // Remove row k from b and row/column k from A.
  Eigen::SparseMatrix<double> newA(n - 1, n - 1);
  Eigen::VectorXd reducedB(n - 1);
  Eigen::VectorXd newX(n - 1);

  for (Eigen::Index i = 0, ii = 0; i < n; ++i) {
    if (i == k)
      continue;

    reducedB[ii] = newB[i];
    newX[ii] = x[i];
    ++ii;
  }

  for (Eigen::Index outer = 0; outer < A.outerSize(); ++outer) {
    for (Eigen::SparseMatrix<double>::InnerIterator it(A, outer); it; ++it) {
      const Eigen::Index r = it.row();
      const Eigen::Index c = it.col();

      if (r == k || c == k)
        continue;

      const Eigen::Index newR = r - (r > k ? 1 : 0);
      const Eigen::Index newC = c - (c > k ? 1 : 0);

      newA.insert(newR, newC) = it.value();
    }
  }

  newA.makeCompressed();

  // Record the fixed value.
  const Eigen::Index fixedIndex = fixedValues.size();
  fixedValues.push_back(value);
  scalarLocation[coordId] = -static_cast<int>(fixedIndex) - 1;

  // Update all remaining mappings.
  for (auto &loc : scalarLocation) {
    if (loc > k)
      --loc;
  }

  A = std::move(newA);
  b = std::move(reducedB);
  x = std::move(newX);
}
