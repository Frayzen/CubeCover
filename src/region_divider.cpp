#include "region_divider.hpp"
#include "constants.hpp"
#include "graph.hpp"
#include "mesh_saver.hpp"
#include <cmath>
#include <iostream>
#include <queue>
#include <set>
#include <string>
RegionDivider::RegionDivider(EdgeGraph &g)
    : edgeGraph(g), mesh(g.mesh), nbRegion(1), regions() {
  for (size_t i = 0; i < g.mesh.nbBoundaryNode; i++)
    regions.set(i, 0);
  for (Axis d = 0; d < DIMENSION; d++)
    loopZones[d].insertZone({0}, 0);
  computeNodeNormals();
}

static Axis getAxis(const Vertex &v) {
  size_t res = 0;
  double maxVal = abs(Vertex(1, 0, 0).dot(v));
  for (size_t d = 1; d < DIMENSION; d++) {
    Vertex a = Vertex::Zero();
    a[d] = 1;
    double cur = abs(v.dot(a));
    if (cur > maxVal) {
      maxVal = cur;
      res = d;
    }
  }
  return res;
}

void RegionDivider::computeNodeNormals() {
  nodeNormals = std::vector<std::array<bool, 3>>(mesh.nbBoundaryNode,
                                                 {false, false, false});
  for (FaceId faceId = 0; faceId < mesh.nbBoundaryFace; faceId++) {
    Face &f = mesh.faces[faceId];
    Vertex normal = f.computeNormal(mesh);
    Axis a = getAxis(normal);
    for (size_t i = 0; i < NB_NODE_PER_FACE; i++)
      nodeNormals[f[i]][a] = true;
  }
}

RegionId RegionDivider::fillConnectedFaces(NodeId startingEdgeId,
                                           std::vector<EdgeId> loopEdgeIds,
                                           size_t axis) {
  size_t targetRegion = regions[mesh.edges[startingEdgeId][0]];
  size_t newRegion = nbRegion;
  nbRegion++;
  std::vector<bool> processed(mesh.nbBoundaryEdge, false);
  std::queue<size_t> toProcess;

  Edge &startingEdge = mesh.edges[startingEdgeId];
  Vertex &firstVertex = mesh.vertices[startingEdge[0]];
  Vertex &secondVertex = mesh.vertices[startingEdge[1]];
  size_t followUpNode = startingEdge[0];
  if (firstVertex[axis] > secondVertex[axis])
    followUpNode = startingEdge[1];

  processed[startingEdgeId] = true;
  for (auto adjEdgeId : edgeGraph.adjacency[startingEdgeId]) {
    if (mesh.edges[adjEdgeId].containsNode(followUpNode)) {
      toProcess.push(adjEdgeId);
      processed[adjEdgeId] = true;
    }
  }
  regions.set(followUpNode, newRegion);

  while (!toProcess.empty()) {
    size_t edgeId = toProcess.front();
    toProcess.pop();
    assert(edgeId < mesh.nbBoundaryEdge);
    Edge &edge = mesh.edges[edgeId];
    assert(regions[edge[0]] == newRegion || regions[edge[1]] == newRegion);
    size_t newEdgeLocal = regions[edge[0]] == newRegion ? 0 : 1;
    size_t followUpNode = edge[newEdgeLocal == 0 ? 1 : 0];
    bool partOfTheLoop = std::find(loopEdgeIds.begin(), loopEdgeIds.end(),
                                   edgeId) != loopEdgeIds.end();
    bool boundaryEdge = regions[followUpNode] != targetRegion;
    if (partOfTheLoop || boundaryEdge)
      followUpNode = edge[newEdgeLocal];
    regions.set(followUpNode, newRegion);
    for (auto adjEdgeId : edgeGraph.adjacency[edgeId]) {
      if (processed[adjEdgeId])
        continue;
      if (mesh.edges[adjEdgeId].containsNode(followUpNode)) {
        toProcess.push(adjEdgeId);
        processed[adjEdgeId] = true;
      }
    }
  }
  return newRegion;
}

void RegionDivider::insertLoop(std::vector<size_t> edgeIds, size_t axis) {
  paths.push_back(std::make_pair(edgeIds, axis));

  Edge &firstEdge = mesh.edges[edgeIds[0]];
  RegionId oldRegion = regions[firstEdge[0]];
  ZoneId oldZone = loopZones[axis].getZone(oldRegion);

  LoopZone loopZone;

  for (size_t i = 0; i < edgeIds.size(); i++) {
    size_t edgeId = edgeIds[i];

    Edge &e = mesh.edges[edgeId];
    RegionId curOldRegion = regions[e[0]];
    bool sameRegion = curOldRegion == regions[e[1]];
    if (!sameRegion)
      continue;
    RegionId newRegion = fillConnectedFaces(edgeId, edgeIds, axis);
    loopZone.push_back(newRegion);
    for (Axis d = 0; d < DIMENSION; d++) {
      if (d == axis)
        continue;
      loopZones[d].addRegionToZone(loopZones[d].getZone(curOldRegion),
                                   newRegion);
    }
    std::cout << "GOT " << nbRegion << " REGIONS " << std::endl;
  }
  loopZones[axis].insertZone(loopZone, oldZone);
  for (size_t i = 0; i < DIMENSION; i++) {
    std::cout << "LOOP ZONE " << i << " GOT " << loopZones[i].size()
              << std::endl;
  }
}

std::vector<NodeId> RegionDivider::getCandidates(RegionId regionId) {
  std::vector<NodeId> candidates;
  ssize_t maxScore = -1;

  // TODO change to retrieve from region
  bool axisNormals[3] = {true, true, true};

  auto nodeIds = regions.getAll(regionId);
  for (NodeId nodeId : nodeIds) {
    ssize_t score = 0;
    for (size_t d = 0; d < DIMENSION; d++) {
      if (axisNormals[d] && nodeNormals[nodeId][d])
        score++;
      if (axisNormals[d] && !nodeNormals[nodeId][d])
        score--;
    }
    if (score > maxScore) {
      candidates.clear();
      maxScore = score;
    }
    if (score == maxScore)
      candidates.push_back(nodeId);
  }
  return candidates;
}

struct CornerCandidate {
  double axisValue;
  NodeId nodeId;
  RegionId regionId;
};

struct Order {
  bool operator()(const CornerCandidate &a, const CornerCandidate &b) const {
    if (a.axisValue == b.axisValue)
      return a.nodeId > b.nodeId;
    return a.axisValue < b.axisValue; // Sorts in ascending order
  }
};

std::vector<double>
AxisLoopZones::computeMids(TetrahedralMesh &mesh,
                           std::vector<std::vector<NodeId>> &regionCandidates,
                           size_t axis) {
  std::vector<double> res(loopZones.size());
  std::vector<size_t> regionToCandidateCount(regionCandidates.size());
  for (size_t zoneId = 0; zoneId < loopZones.size(); zoneId++) {
    const LoopZone &zone = loopZones[zoneId];
    std::set<CornerCandidate, Order> zoneCandidates;
    for (size_t i = 0; i < zone.size(); i++) {
      RegionId regionId = zone[i];
      for (NodeId candidateId : regionCandidates[regionId]) {
        CornerCandidate c{
            .axisValue = mesh.vertices[candidateId][axis],
            .nodeId = candidateId,
            .regionId = regionId,
        };
        zoneCandidates.insert(c);
      }
      regionToCandidateCount[regionId] = regionCandidates[regionId].size();
    }
    assert(!zoneCandidates.empty());
    auto leftPtr = zoneCandidates.begin();
    auto rightPtr = zoneCandidates.end();
    rightPtr--;
    while (rightPtr != leftPtr) {
      const CornerCandidate &cur = *rightPtr;
      if (regionToCandidateCount[cur.regionId] == 1)
        break;
      regionToCandidateCount[cur.regionId]--;
      rightPtr--;
    }
    while (leftPtr != rightPtr) {
      const CornerCandidate &cur = *leftPtr;
      if (regionToCandidateCount[cur.regionId] == 1)
        break;
      regionToCandidateCount[cur.regionId]--;
      leftPtr++;
    }
    double mid = (leftPtr->axisValue + rightPtr->axisValue) / 2.0;
    res[zoneId] = mid;
  }
  return res;
}

std::vector<NodeId> RegionDivider::computeCorners() {
  std::vector<std::vector<NodeId>> regionCandidates(nbRegion);
  for (RegionId regionId = 0; regionId < nbRegion; regionId++)
    regionCandidates[regionId] = getCandidates(regionId);
  std::array<std::vector<double>, DIMENSION> zoneMidPerAxis;
  for (Axis d = 0; d < DIMENSION; d++)
    zoneMidPerAxis[d] = loopZones[d].computeMids(mesh, regionCandidates, d);

  std::vector<NodeId> regionCorners(nbRegion);
  for (RegionId regionId = 0; regionId < nbRegion; regionId++) {
    Vertex target;
    for (Axis d = 0; d < DIMENSION; d++) {
      ZoneId zoneId = loopZones[d].getZone(regionId);
      target[d] = zoneMidPerAxis[d][zoneId];
    }
    double closestDist = INFINITY;
    NodeId closestId = 0;
    for (NodeId nodeId : regionCandidates[regionId]) {
      double distance = (mesh.vertices[nodeId] - target).norm();
      if (distance < closestDist) {
        closestId = nodeId;
        closestDist = distance;
      }
    }
    regionCorners[regionId] = closestId;
  }
  return regionCorners;
}

void RegionDivider::linkCorners(NodeId start, NodeId end) {
  Graph g(mesh.nbBoundaryNode);
  AdjacencyWeight weights(mesh.nbBoundaryNode);
  for (EdgeId edgeId = 0; edgeId < mesh.nbBoundaryEdge; edgeId++) {
    Edge &e = mesh.edges[edgeId];
    g.adjacency[e[0]].push_back(e[1]);
    g.adjacency[e[0]].push_back(e[1]);
  }
}

void RegionDivider::computePaths(std::vector<NodeId> corners) {
  assert(corners.size() == nbRegion);
  for (RegionId startRegionId = 0; startRegionId < nbRegion; startRegionId++) {
    std::vector<size_t> sharedLoops(nbRegion, 0);
    for (Axis d = 0; d < DIMENSION; d++) {
      ZoneId zoneId = loopZones[d].getZone(startRegionId);
      for (RegionId other : loopZones[d][zoneId])
        sharedLoops[other]++;
    }
    for (RegionId endRegionId = 0; endRegionId < nbRegion; endRegionId++) {
      if (endRegionId == startRegionId)
        continue;
      assert(sharedLoops[endRegionId] < 3);
      if (sharedLoops[endRegionId] < 2 || endRegionId > startRegionId)
        continue;
      linkCorners(corners[startRegionId], corners[endRegionId]);
    }
  }
}

void RegionDivider::saveMesh(std::string filename) {
  std::vector<NodeId> corners = computeCorners();
  std::vector<bool> isCorner(mesh.vertices.size(), false);

  std::vector<ZoneId> xZone(mesh.vertices.size(), false);
  std::vector<ZoneId> yZone(mesh.vertices.size(), false);
  std::vector<ZoneId> zZone(mesh.vertices.size(), false);
  for (NodeId nodeId = 0; nodeId < mesh.nbBoundaryNode; nodeId++) {
    RegionId regionId = regions[nodeId];
    xZone[nodeId] = loopZones[0].getZone(regionId);
    yZone[nodeId] = loopZones[1].getZone(regionId);
    zZone[nodeId] = loopZones[2].getZone(regionId);
  }

  for (NodeId nodeId : corners)
    isCorner[nodeId] = true;
  std::vector<int> colors(mesh.vertices.size(), -1);
  for (size_t nodeId = 0; nodeId < mesh.nbBoundaryNode; nodeId++)
    colors[nodeId] = static_cast<size_t>(getFaceRegion(nodeId));

  MeshSaver saver(mesh.vertices, mesh.edges, mesh.tets, filename);
  saver.addField<ElementType::Node>("region", colors);
  saver.addField<ElementType::Node>("isCorner", isCorner);
  saver.addField<ElementType::Node>("xZone", xZone);
  saver.addField<ElementType::Node>("yZone", yZone);
  saver.addField<ElementType::Node>("zZone", zZone);
  saver.saveMesh();
}
