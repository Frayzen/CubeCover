#include "mesh_loader.hpp"
#include "signed_axis.hpp"
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
  for (Axis d = 0; d < DIMENSION; d++) {
    size_t nbSlice = mesh.sliceManagers[d].nbSlice();
    mesh.slice(d, nbSlice / 5);
    mesh.slice(d, 4 * nbSlice / 5);
  }
  mesh.splitByConnectedRegions();

  VoxelMeshSaver saver(mesh);

  std::vector<size_t> colors(mesh.voxels.size());
  std::vector<double> xNormal(mesh.voxels.size());
  std::vector<double> yNormal(mesh.voxels.size());
  std::vector<double> zNormal(mesh.voxels.size());
  size_t regionId = 0;
  for (Region region : mesh.regions) {
    for (VoxelId id : region) {
      colors[id] = regionId;
      xNormal[id] = getIntegerSignValue(region.signedNormals[0]);
      yNormal[id] = getIntegerSignValue(region.signedNormals[1]);
      zNormal[id] = getIntegerSignValue(region.signedNormals[2]);
    }
    regionId++;
  }
  saver.addField("region", colors);
  saver.addField("xNormal", xNormal);
  saver.addField("yNormal", yNormal);
  saver.addField("zNormal", zNormal);

  saver.saveMesh();
}
