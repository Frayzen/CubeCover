#include "constants.hpp"
#include "frame_field.hpp"
#include "iterator.hpp"
#include "mesh_loader.hpp"
#include "mesh_saver.hpp"
#include "span_tree.hpp"
#include "trilaterate.hpp"
#include <Eigen/Core>
#include <Eigen/IterativeLinearSolvers>

#include <Eigen/Eigen>
#include <Eigen/Eigenvalues>
#include <iostream>

struct BBox {
  Vertex max;
  Vertex min;
};

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
  std::vector<Eigen::Matrix3d> frames(mesh.tets.size());
  fillRandomPermutations(frames, mesh);

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

  CubeCoverIterator system(mesh, frames);
  system.initialize();
  system.solve();

  // size_t dimensionOffset = NB_UNKNOWN_PER_TET;
  // auto setLocalVertices = [&](size_t tetId, size_t localNodeId, Vertex value)
  // {
  //   size_t offset = tetId * NB_UNKNOWN_PER_TET;
  //   for (size_t dim = 0; dim < DIMENSION; dim++)
  //     system.x[offset + dim * NB_NODE_PER_TET + localNodeId] = value[dim];
  // };
  // auto getLocalVertices = [&](size_t tetId, size_t localNodeId) -> Vertex {
  //   size_t offset = tetId * NB_UNKNOWN_PER_TET;
  //   Vertex res;
  //   for (size_t dim = 0; dim < DIMENSION; dim++)
  //     res[dim] = system.x[offset + dim * NB_NODE_PER_TET + localNodeId];
  //   return res;
  // };
  // auto getGlobalVertices = [&](size_t tetId, size_t globalNodeId) {
  //   auto tetNodes = mesh.tets[tetId];
  //   for (size_t i = 0; i < tetNodes.size(); i++) {
  //     if (tetNodes[i] == globalNodeId)
  //       return getLocalVertices(tetId, i);
  //   }
  //   assert(false);
  // };
  // auto setGlobalVertices = [&](size_t tetId, size_t globalNodeId,
  //                              Vertex value) {
  //   auto tetNodes = mesh.tets[tetId];
  //   for (size_t i = 0; i < tetNodes.size(); i++) {
  //     if (tetNodes[i] == globalNodeId)
  //       return setLocalVertices(tetId, i, value);
  //   }
  //   assert(false);
  // };

  // std::vector<Eigen::Vector3d> vertices(mesh.tets.size() * NB_NODE_PER_TET);
  // for (size_t i = 0; i < mesh.tets.size(); i++) {
  //   for (size_t j = 0; j < NB_NODE_PER_TET; j++)
  //     vertices[NB_NODE_PER_TET * i + j] = getLocalVertices(i, j);
  // }
  // // saveMeshVtk(vertices, "nospanning.vtk");
  // std::cout << "SAVED NO SPANNING" << std::endl;

  // std::vector<bool> include(mesh.tets.size());
  // std::fill(include.begin(), include.end(), false);
  // for (size_t faceId : spanTree) {
  //   Face &face = mesh.faces[faceId];
  //   auto faceNodes = face.vertices;
  //   if (face.adjacentCells[1] == UNDEFINED_ELEMENT)
  //     continue;
  //   size_t parent = face.adjacentCells[0];
  //   size_t child = face.adjacentCells[1];
  //   include[parent] = true;
  //   include[child] = true;

  //   // Find the child tip
  //   auto childNodes = mesh.tets[child];
  //   size_t childTipLocalId = NB_NODE_PER_TET;
  //   for (childTipLocalId = 0; childTipLocalId < NB_NODE_PER_TET;
  //        childTipLocalId++) {
  //     size_t tipGlobalId = childNodes[childTipLocalId];
  //     if (std::find(faceNodes.begin(), faceNodes.end(), tipGlobalId) ==
  //         faceNodes.end())
  //       break;
  //   }
  //   assert(childTipLocalId < NB_NODE_PER_TET);

  //   size_t parentTipLocalId = NB_NODE_PER_TET;
  //   auto parentNodes = mesh.tets[parent];
  //   for (parentTipLocalId = 0; parentTipLocalId < NB_NODE_PER_TET;
  //        parentTipLocalId++) {
  //     size_t tipGlobalId = parentNodes[parentTipLocalId];
  //     if (std::find(faceNodes.begin(), faceNodes.end(), tipGlobalId) ==
  //         faceNodes.end())
  //       break;
  //   }
  //   assert(parentTipLocalId < NB_NODE_PER_TET);

  //   Vertex childTip(getLocalVertices(child, childTipLocalId));

  //   Vertex oldFace[NB_NODE_PER_FACE];
  //   Vertex newFace[NB_NODE_PER_FACE];

  //   for (size_t i = 0; i < NB_NODE_PER_FACE; i++) {
  //     size_t node = faceNodes[i];
  //     Vertex fp = getGlobalVertices(parent, node);
  //     Vertex fc = getGlobalVertices(child, node);
  //     oldFace[i] = fc;
  //     newFace[i] = fp;
  //     setGlobalVertices(child, node, fp);
  //   }

  //   Vertex &parentTip = mesh.vertices[parentNodes[parentTipLocalId]];
  //   auto newTip = tipProject(oldFace, newFace, childTip);
  //   setLocalVertices(child, childTipLocalId, newTip);
  // }

  // vertices.clear();
  // for (size_t i = 0; i < mesh.tets.size(); i++) {
  //   if (!include[i])
  //     continue;
  //   for (size_t j = 0; j < NB_NODE_PER_TET; j++)
  //     vertices.push_back(getLocalVertices(i, j));
  // }
  // saveMeshVtk(vertices, frames);

  // vertices.clear();
  // for (size_t i = 0; i < mesh.tets.size(); i++) {
  //   for (size_t j = 0; j < NB_NODE_PER_TET; j++)
  //     vertices.push_back(mesh.vertices[mesh.tets[i][j]]);
  // }

  // saveMeshVtk(vertices, frames, "input.vtk");

  return EXIT_SUCCESS;
}
