#include "dual_saver.hpp"
#include "ids.hpp"
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

  std::vector<Voxel> junctions;
// #define COMPUTEJUNCTIONS
#ifdef COMPUTEJUNCTIONS
  // #include "junctions.hpp"
  SurfacePatches patches = buildSurfacePatches(tetMesh);
  std::vector<size_t> patchIds(tetMesh.nbBoundaryFace);
  std::fill(patchIds.begin(), patchIds.end(), 0);
  size_t curPatch = 0;
  for (auto patch : patches.patches) {
    for (FaceId f : patch) {
      patchIds[f] = curPatch;
      assert(f < tetMesh.nbBoundaryFace);
    }
    curPatch++;
  }

  junctions = computeMedialAxisJunctions(mesh, patches);
  std::cout << "FOUND " << junctions.size() << " JUNCTIONS " << std::endl;

  MeshSaver tetMeshSaver(tetMesh, "patches.vtu");
  tetMeshSaver.addField<ElementType::FaceBoundaryType>("patchIds", patchIds);
  tetMeshSaver.saveMesh();
  tetMeshSaver.saveFaces();
#endif

  Voxel mid = mesh.bounds() / 2;
  for (Axis d = 0; d < DIMENSION; d++) {
    mesh.cut(d, mid);
  }

  mesh.cut(0, mesh.bounds() / 7);
  mesh.cut(0, (6 * mesh.bounds()) / 7);
  mesh.cut(2, mesh.bounds() / 4);

  Voxel cut = (3 * mesh.bounds()) / 4;
  cut.z() -= 11;
  std::cout << "HERE " << std::endl;
  mesh.cut(2, cut);

  for (Axis d = 0; d < DIMENSION; d++) {
    mesh.saveDotfile(d, "region_adjacency" + std::to_string(d) + ".dot");
  }

  std::vector<Voxel> corners;
  std::vector<size_t> pathIds(mesh.nbVoxel);
  std::fill(pathIds.begin(), pathIds.end(), 0);

  saveDual("dual.vtu", mesh);
// #define COMPUTEPATHS
#ifdef COMPUTEPATHS
  corners = mesh.findCorners();
  std::vector<std::vector<Voxel>> paths = mesh.findPaths();

  int pathId = 1;
  for (auto path : paths) {
    for (Voxel v : path)
      pathIds[mesh.getId(v)] = pathId;
    pathId++;
  }
#endif

  VoxelMeshSaver saver(mesh);

  std::vector<bool> junctionFlags(mesh.nbVoxel);
  std::fill(junctionFlags.begin(), junctionFlags.end(), false);
  for (auto voxel : junctions)
    junctionFlags[mesh.getId(voxel)] = true;

  std::vector<bool> cornersFlag(mesh.nbVoxel);

  std::vector<bool> allowedCutX(mesh.nbVoxel);
  std::vector<bool> allowedCutY(mesh.nbVoxel);
  std::vector<bool> allowedCutZ(mesh.nbVoxel);

  std::vector<bool> boundaryFlag(mesh.nbVoxel);

  std::fill(cornersFlag.begin(), cornersFlag.end(), false);
  std::vector<size_t> colors(mesh.nbVoxel);
  std::fill(colors.begin(), colors.end(), NO_REGION);

  std::vector<double> xNormal(mesh.nbVoxel);
  std::vector<double> yNormal(mesh.nbVoxel);
  std::vector<double> zNormal(mesh.nbVoxel);

  for (size_t i = 0; i < mesh.voxels.dimension(0); i++) {
    for (size_t j = 0; j < mesh.voxels.dimension(1); j++) {
      for (size_t k = 0; k < mesh.voxels.dimension(2); k++) {
        Voxel v(i, j, k);
        VoxelId id = mesh.getId(v);
        if (id == NO_VOXEL)
          continue;
        RegionId regionId = mesh.voxelToRegion(i, j, k);
        Region &region = mesh.regions[regionId];
        colors[id] = regionId;
        xNormal[id] = getIntegerSignValue(region.signedNormals[0]);
        yNormal[id] = getIntegerSignValue(region.signedNormals[1]);
        zNormal[id] = getIntegerSignValue(region.signedNormals[2]);
        boundaryFlag[id] = mesh.voxelBoundaryFlag(i, j, k);

        allowedCutX[id] = mesh.isCutAllowed(0, v);
        allowedCutY[id] = mesh.isCutAllowed(1, v);
        allowedCutZ[id] = mesh.isCutAllowed(2, v);
#ifdef COMPUTEPATHS
        if (v == corners[regionId])
          cornersFlag[id] = true;
#endif
      }
    }
  }

  saver.addField("boundaryFlag", boundaryFlag);
  saver.addField("region", colors);
  saver.addField("xNormal", xNormal);
  saver.addField("yNormal", yNormal);
  saver.addField("zNormal", zNormal);
  saver.addField("corner", cornersFlag);
  saver.addField("pathIds", pathIds);
  saver.addField("junctionFlag", junctionFlags);
  saver.addField("allowedCutX", allowedCutX);
  saver.addField("allowedCutY", allowedCutY);
  saver.addField("allowedCutZ", allowedCutZ);

  saver.saveMesh();
}
