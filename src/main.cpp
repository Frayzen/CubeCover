#include "constants.hpp"
#include "frame_field.hpp"
#include "iterator.hpp"
#include "mesh_loader.hpp"
#include "mesh_saver.hpp"
#include <Eigen/Core>
#include <Eigen/IterativeLinearSolvers>

#include <Eigen/Eigen>
#include <Eigen/Eigenvalues>
#include <exception>
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
  // for (size_t i = 0; i < mesh.vertices.size(); i++)
  //   mesh.vertices[i] *= 100;

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

  std::vector<Vertex> vertices;
  std::vector<Eigen::Vector3i> fixedMask;
  for (size_t i = 0; i < mesh.tets.size(); i++) {
    for (size_t j = 0; j < NB_NODE_PER_TET; j++) {
      vertices.push_back(mesh.vertices[mesh.tets[i][j]]);
      fixedMask.push_back(Eigen::Vector3i::Zero());
    }
  }
  saveMeshVtk(vertices, frames, fixedMask, "input.vtu");

  CubeCoverIterator system(mesh, frames);
  system.initialize();

  size_t totalIter = system.coordsToRemove.size();
  std::vector<std::string> meshFiles;
  size_t i = 0;

  while (!system.coordsToRemove.empty()) {
    std::cout << "============ ITERATION " << i << " / " << totalIter
              << " ============" << std::endl;

    try {
      system.iterate();
    } catch (std::exception &e) {
      std::cout << "Error " << e.what() << std::endl;
      saveMeshCollection("output.pvd", meshFiles);
      exit(1);
    }

    std::string filename = "output_" + std::to_string(i) + ".vtu";
    system.saveMesh(filename);
    meshFiles.push_back(filename);
    saveMeshCollection("output.pvd", meshFiles);

    ++i;
    exit(1);
  }

  saveMeshCollection("output.pvd", meshFiles);

  return EXIT_SUCCESS;
}
