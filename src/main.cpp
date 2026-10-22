#include "constants.hpp"
#include "frame_field.hpp"
#include "iterator.hpp"
#include "mesh_loader.hpp"
#include "mesh_saver.hpp"
#include "span_tree.hpp"
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

  system.applySpanTree(spanTree);
  system.saveMesh("output.vtk");

  std::vector<Vertex> vertices;
  for (size_t i = 0; i < mesh.tets.size(); i++) {
    for (size_t j = 0; j < NB_NODE_PER_TET; j++)
      vertices.push_back(mesh.vertices[mesh.tets[i][j]]);
  }
  saveMeshVtk(vertices, frames, "input.vtk");

  return EXIT_SUCCESS;
}
