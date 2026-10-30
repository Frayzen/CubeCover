#ifndef AXIS_GRAPH_HPP
#define AXIS_GRAPH_HPP

#include "edge_graph.hpp"
#include "many_to_one.hpp"
#include <stdexcept>

using NodeId = size_t;
using FaceId = size_t;
using EdgeId = size_t;
using RegionId = size_t;
using Axis = size_t;
using Path = std::pair<std::vector<EdgeId>, Axis>;
using ZoneId = size_t;

class RegionGraph;

using LoopZone = std::vector<RegionId>;
class AxisLoopZones {
  std::vector<LoopZone> loopZones;

public:
  std::vector<double>
  computeMids(TetrahedralMesh &mesh,
              std::vector<std::vector<NodeId>> &regionCandidates, size_t axis);
  const LoopZone &operator[](ZoneId zoneId) const { return loopZones[zoneId]; };

  ZoneId getZone(RegionId region) const {
    for (ZoneId i = 0; i < loopZones.size(); i++) {
      auto r = loopZones[i];
      if (std::find(r.begin(), r.end(), region) != r.end())
        return i;
    }
    throw std::invalid_argument("Invalid region");
  };

  void addRegionToZone(ZoneId zoneId, RegionId regionId) {
    loopZones[zoneId].push_back(regionId);
  }
  size_t size() { return loopZones.size(); }

  void insertZone(LoopZone loopZone, ZoneId position) {
    loopZones.insert(loopZones.begin() + position, loopZone);
  }
};

class RegionGraph {
private:
  EdgeGraph &g;
  TetrahedralMesh &mesh;
  size_t nbRegion;
  ManyToOne<NodeId, RegionId> regions;
  AxisLoopZones loopZones[DIMENSION];

  std::vector<std::array<bool, DIMENSION>> nodeNormals;
  void computeNodeNormals();

  RegionId fillConnectedFaces(size_t startingFace, std::vector<size_t> edgeIds,
                              size_t axis);

public:
  std::vector<Path> paths;
  RegionGraph(EdgeGraph &g);

  void insertLoop(std::vector<size_t> edgesId, size_t axis);
  size_t getFaceRegion(size_t nodeId) { return regions[nodeId]; };
  std::vector<NodeId> computeCorners();

  std::vector<NodeId> getCandidates(RegionId regionId);

  void saveMesh(std::string filename);
  friend AxisLoopZones;
};

#endif
