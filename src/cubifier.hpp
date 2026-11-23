#ifndef CUBIFIER_HPP
#define CUBIFIER_HPP

#include "constants.hpp"
#include "ids.hpp"
#include "surface_mesh.hpp"

using Hex = std::array<NodeId, NB_NODE_PER_HEX>;

class Cubifier {
private:
  void build();
  void applyConstraints(const Eigen::Tensor<NodeId, DIMENSION> &nodeIds,
                        const Voxel &start, std::vector<bool> &constrainedFlag,
                        bool last);

public:
  std::vector<Hex> hexes;
  std::vector<Vertex> vertices;
  std::vector<bool> constrainedFlag;

  const SurfaceMesh &surfaceMesh;
  const VoxelMesh &mesh;
  const Region &region;

  std::vector<Voxel> nodeGridCoordinates;

  Cubifier(const SurfaceMesh &surfaceMesh, RegionId regionId)
      : surfaceMesh(surfaceMesh), mesh(surfaceMesh.m),
        region(mesh.regions[regionId]) {
    build();
  }

  void exportMesh(std::string filename) const;
};

#endif
