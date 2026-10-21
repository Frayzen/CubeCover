#include "process.hpp"
#include "constants.hpp"
#include "mesh_saver.hpp"
#include "trilaterate.hpp"
#include "vertex.hpp"
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SparseCore>
#include <iostream>

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

static void removeSystemRowAndCol(CubeCoverParameters &params,
                                  size_t idToRemove) {
  Eigen::SparseMatrix<double> &A = params.system;
  Eigen::VectorXd &b = params.knowns;
  Eigen::VectorXd &x = params.result;

  const Eigen::Index n = A.rows();
  const Eigen::Index j = static_cast<Eigen::Index>(idToRemove);

  if (j >= n) {
    throw std::out_of_range("Index to remove is out of bounds");
  }

  double removedX = x[idToRemove];

  // Rebuild the matrix without row j and column j.
  Eigen::SparseMatrix<double> reducedA(n - 1, n - 1);
  std::vector<Eigen::Triplet<double>> triplets;
  Eigen::VectorXd reducedX(n - 1);
  Eigen::VectorXd reducedB(n - 1);
  triplets.reserve(static_cast<size_t>(A.nonZeros()));

  for (Eigen::Index col = 0; col < A.outerSize(); ++col) {
    if (col == j)
      continue;
    const Eigen::Index newCol = col < j ? col : col - 1;
    for (Eigen::SparseMatrix<double>::InnerIterator it(A, col); it; ++it) {
      if (it.row() == j)
        continue;
      const Eigen::Index row = it.row() < j ? it.row() : it.row() - 1;
      triplets.emplace_back(row, newCol, it.value());
      reducedB[newCol] -= it.value() * removedX;
      reducedX[newCol] = x[col];
    }
  }
  reducedA.setFromTriplets(triplets.begin(), triplets.end());
  A.swap(reducedA);
  b.swap(reducedB);
  x.swap(reducedX);
}

void cubeCoverIteration(CubeCoverParameters &params, size_t rowToRemove) {

  Eigen::ConjugateGradient<Eigen::SparseMatrix<double>,
                           Eigen::Lower | Eigen::Upper>
      solver;
  TetrahedralMesh &mesh = params.mesh;
  solver.setMaxIterations(10000);
  solver.setTolerance(1e-3);
  solver.compute(params.system);

  if (solver.info() != Eigen::Success) {
    std::cerr << "Failed to initialize CG solver\n";
  } else {
    std::cout << "Begin solving" << std::endl;
    params.unspannedSolution =
        solver.solveWithGuess(params.knowns, params.result);
    std::cout << "End solving" << std::endl;

    std::cerr << "Iterations: " << solver.iterations() << '\n';
    std::cerr << "Estimated error: " << solver.error() << '\n';

    if (solver.info() != Eigen::Success) {
      std::cerr << "CG failed to converge\n";
    }
  }

  params.result = params.unspannedSolution;
  size_t dimensionOffset = NB_UNKNOWN_PER_TET;
  auto setLocalVertices = [&](size_t tetId, size_t localNodeId, Vertex value) {
    size_t offset = tetId * NB_UNKNOWN_PER_TET;
    for (size_t dim = 0; dim < DIMENSION; dim++)
      params.result[offset + dim * NB_NODE_PER_TET + localNodeId] = value[dim];
  };
  auto getLocalVertices = [&](size_t tetId, size_t localNodeId) -> Vertex {
    size_t offset = tetId * NB_UNKNOWN_PER_TET;
    Vertex res;
    for (size_t dim = 0; dim < DIMENSION; dim++)
      res[dim] = params.result[offset + dim * NB_NODE_PER_TET + localNodeId];
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

  // for (size_t i = 0; i < mesh.tets.size(); i++) {
  //   for (size_t j = 0; j < NB_NODE_PER_TET; j++)
  //     vertices[NB_NODE_PER_TET * i + j] = getLocalVertices(i, j);
  // }
  // saveMeshVtk(vertices, "nospanning.vtk");
  // std::cout << "SAVED NO SPANNING" << std::endl;

  for (size_t faceId : params.spanTree) {
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

    Vertex &parentTip = mesh.vertices[parentNodes[parentTipLocalId]];
    auto newTip = tipProject(oldFace, newFace, childTip);
    setLocalVertices(child, childTipLocalId, newTip);
  }

  // Face &boundaryFace = mesh.faces[boundaryFaceId];
  // size_t tetId = boundaryFace.adjacentCells[0];
  // auto tet = mesh.tets[tetId];
  // assert(boundaryFace.adjacentCells[1] == UNDEFINED_ELEMENT);

  // size_t offset = tetId * NB_UNKNOWN_PER_TET;

  // Vertex textureNodes[NB_NODE_PER_FACE];
  // for (size_t i = 0; i < NB_NODE_PER_FACE; i++)
  //   textureNodes[i] = getGlobalVertices(tetId, boundaryFace.vertices[i]);
  // Vertex normal = getNormal(textureNodes);
  // size_t axis = getPrincipalAxis(normal);

  // std::vector<size_t> rowsToRemove;

  // for (size_t nodeId : boundaryFace.vertices) {
  //   size_t tetLocalId = 0;
  //   while (tetLocalId < NB_NODE_PER_TET && tet[tetLocalId] != nodeId)
  //     tetLocalId++;
  //   assert(tetLocalId < NB_NODE_PER_TET);
  //   for (size_t i = 0; i < DIMENSION; i++) {
  //     size_t resultId = offset + i * NB_NODE_PER_TET + tetLocalId;
  //     double previousValue = params.result[resultId];
  //     params.result[resultId] = std::round(previousValue);
  //     rowsToRemove.push_back(resultId);
  //   }
  // }
  std::cout << "REACH REMOVE" << std::endl;
  removeSystemRowAndCol(params, rowToRemove);
}
