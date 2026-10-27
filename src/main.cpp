#include "constants.hpp"
#include "iterator.hpp"
#include "mesh_loader.hpp"
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

  Eigen::Matrix3d frame;
  frame.col(0) << 1, 0, -0.2;
  frame.col(1) << 0, 1, -0.2;
  frame.col(2) << 0, 0, 1;

  std::vector<Eigen::Matrix3d> frames(mesh.tets.size(), frame);

  // fillRandomPermutations(frames, mesh);

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

  std::vector<std::string> meshFiles;

  try {
    system.iterate();
  } catch (std::exception &e) {
    std::cout << "Error " << e.what() << std::endl;
    exit(1);
  }

  std::string filename = "output.vtu";
  system.saveMesh(filename);

  return EXIT_SUCCESS;
}
