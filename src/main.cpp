#include "dual_saver.hpp"
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
#include "junctions.hpp"
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

  std::vector<Voxel> junctions = computeMedialAxisJunctions(mesh, patches);
  std::cout << "FOUND " << junctions.size() << " JUNCTIONS " << std::endl;

  MeshSaver tetMeshSaver(tetMesh, "patches.vtu");
  tetMeshSaver.addField<ElementType::FaceBoundaryType>("patchIds", patchIds);
  tetMeshSaver.saveMesh();
  tetMeshSaver.saveFaces();
#endif

  for (Axis d = 0; d < DIMENSION; d++) {
    double nbSlice = mesh.voxels.dimension(d);
    mesh.cut(d, nbSlice / 2);
  }

  saveDual("dual.vtu", mesh);

  for (Axis d = 0; d < DIMENSION; d++)
    mesh.regionGraphs[d].exportDot("graph_dim" + std::to_string(d) + ".dot");
  std::vector<int> pathIds(mesh.nbVoxel);

  std::vector<Voxel> corners = mesh.findCorners();
  std::vector<std::vector<Voxel>> paths = mesh.findPaths();

  std::fill(pathIds.begin(), pathIds.end(), 0);
  int pathId = 1;
  for (auto path : paths) {
    for (Voxel v : path)
      pathIds[mesh.getId(v)] = pathId;
    pathId++;
  }

  VoxelMeshSaver saver(mesh);

  std::vector<bool> junctionFlags(mesh.nbVoxel);
  std::fill(junctionFlags.begin(), junctionFlags.end(), false);
  for (auto voxel : junctions)
    junctionFlags[mesh.getId(voxel)] = true;

  std::vector<bool> cornersFlag(mesh.nbVoxel);
  std::fill(cornersFlag.begin(), cornersFlag.end(), false);
  std::vector<size_t> colors(mesh.nbVoxel);
  std::fill(colors.begin(), colors.end(), NO_REGION);

  std::vector<double> xNormal(mesh.nbVoxel);
  std::vector<double> yNormal(mesh.nbVoxel);
  std::vector<double> zNormal(mesh.nbVoxel);

  std::vector<double> xLayer(mesh.nbVoxel);
  std::vector<double> yLayer(mesh.nbVoxel);
  std::vector<double> zLayer(mesh.nbVoxel);
  RegionId regionId = 0;
  for (Region region : mesh.regions) {
    for (Voxel v : region) {
      VoxelId id = mesh.getId(v);
      colors[id] = regionId;
      xNormal[id] = getIntegerSignValue(region.signedNormals[0]);
      yNormal[id] = getIntegerSignValue(region.signedNormals[1]);
      zNormal[id] = getIntegerSignValue(region.signedNormals[2]);

      xLayer[id] = region.layerIds[0];
      yLayer[id] = region.layerIds[1];
      zLayer[id] = region.layerIds[2];
      if (v == corners[regionId])
        cornersFlag[id] = true;
    }
    regionId++;
  }
  saver.addField("region", colors);
  saver.addField("xNormal", xNormal);
  saver.addField("yNormal", yNormal);
  saver.addField("zNormal", zNormal);
  saver.addField("xLayer", xLayer);
  saver.addField("yLayer", yLayer);
  saver.addField("zLayer", zLayer);
  saver.addField("corner", cornersFlag);
  saver.addField("pathIds", pathIds);
  saver.addField("junctionFlag", junctionFlags);

  saver.saveMesh();
}
