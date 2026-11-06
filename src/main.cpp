#include "mesh_loader.hpp"
#include "voxel_mesh.hpp"
#include "voxel_mesh_saver.hpp"
#include <Eigen/Core>
#include <Eigen/IterativeLinearSolvers>

#include <Eigen/Eigen>
#include <Eigen/Eigenvalues>
#include <cmath>
#include <cstdlib>
#include <iostream>

int main(int argc, char *argv[]) {
  if (argc != 2) {
    std::cerr << "Usage: " << argv[0] << " <mesh-file.vtu|mesh-file.vtk>\n";
    return EXIT_FAILURE;
  }

  TetrahedralMesh tetMesh;
  if (!loadTetrahedralMesh(argv[1], tetMesh)) {
    return EXIT_FAILURE;
  }

  VoxelMesh mesh(tetMesh);
  for (Axis d = 0; d < DIMENSION; d++)
    mesh.split(d, mesh.sliceManagers[d].nbSlice() / 2);

  VoxelMeshSaver saver(mesh);

  std::vector<size_t> colors(mesh.voxels.size());
  size_t regionId = 0;
  for (Region region : mesh.regions) {
    for (VoxelId id : region)
      colors[id] = regionId;
    regionId++;
  }
  saver.addField("region", colors);

  saver.saveMesh();
}
