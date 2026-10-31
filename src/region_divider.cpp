#include "region_divider.hpp"
#include "constants.hpp"
#include "graph.hpp"
#include "mesh_saver.hpp"
#include "signed_axis.hpp"
#include <Eigen/Core>
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
  regionToNormals.push_back({UNDEFINED, UNDEFINED, UNDEFINED});
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

void RegionDivider::addRegionNormals(RegionId targetRegion, Axis axis) {
  RegionNormals &targetNormals = regionToNormals[targetRegion];
  RegionNormals res = {UNDEFINED, UNDEFINED, UNDEFINED};
  for (Axis d = 0; d < DIMENSION; d++) {
    if (d == axis) {
      if (targetNormals[d] == UNDEFINED) {
        targetNormals[d] = POSITIVE;
        res[d] = NEGATIVE;
      } else if (targetNormals[d] == NEGATIVE) {
        res[d] = NEGATIVE;
        targetNormals[d] = ZERO;
      } else if (targetNormals[d] == POSITIVE)
        res[d] = ZERO;
      continue;
    }
    res[d] = targetNormals[d];
  }
  regionToNormals.push_back(res);
}

RegionId RegionDivider::fillConnectedFaces(NodeId startingEdgeId,
                                           std::vector<EdgeId> loopEdgeIds,
                                           Axis axis) {
  RegionId targetRegion = regions[mesh.edges[startingEdgeId][0]];
  addRegionNormals(targetRegion, axis);
  size_t newRegion = nbRegion;
  nbRegion++;
  std::vector<bool> processed(mesh.nbBoundaryEdge, false);
  std::queue<size_t> toProcess;

  // newRegion is the smallest value
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
      // TODO take into account the sign of the axis normal
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

std::vector<size_t> RegionDivider::linkCorners(NodeId start, NodeId end,
                                               Axis a) {
  Graph graph(mesh.nbBoundaryNode);
  AdjacencyWeight weights(mesh.nbBoundaryNode);
  RegionId regionId = regions[start];
  Axis firstAxis = (a + 1) % 3;
  Axis secondAxis = (a + 2) % 3;
  // bool increase = mesh.vertices[end][a] > mesh.vertices[start][a];
  // if (firstAxis > secondAxis)
  //   std::swap(firstAxis, secondAxis);
  Vertex axis[2] = {
      // left one
      getSignedAxisVertex(regionToNormals[regionId][firstAxis], firstAxis),
      // right one
      getSignedAxisVertex(regionToNormals[regionId][secondAxis], secondAxis),
  };
  for (EdgeId edgeId = 0; edgeId < mesh.nbBoundaryEdge; edgeId++) {
    Edge &e = mesh.edges[edgeId];
    std::array<ssize_t, 2> faces = mesh.edgeToBoundaryFaces[edgeId];
    std::array<Vertex, 2> normals = {
        mesh.faces[faces[0]].computeNormal(mesh),
        mesh.faces[faces[1]].computeNormal(mesh),
    };
    double distance = e.computeLength(mesh);
    NodeId firstNode = e[0];
    NodeId secondNode = e[1];
    // bool firstHigher =
    //     mesh.vertices[firstNode][a] > mesh.vertices[secondNode][a];
    // if (increase ^ firstHigher) // xor
    //   std::swap(firstNode, secondNode);
    graph.adjacency[firstNode].push_back(secondNode);
    graph.adjacency[secondNode].push_back(firstNode);
    double lowestPenality = INFINITY;
    for (size_t i = 0; i < 2; i++) {
      double angle1 = std::acos(normals[i].dot(axis[0]));
      double penality1 =
          angle1 <= 0.955 ? 0 : std::pow((1 + angle1 - 0.955), 2);
      double angle2 = std::acos(normals[(i + 1) % 2].dot(axis[1]));
      double penality2 =
          angle2 <= 0.955 ? 0 : std::pow((1 + angle2 - 0.955), 2);
      double penality = (penality1 + penality2) / 2;
      lowestPenality = std::min(lowestPenality, penality);
    }
    weights[firstNode].push_back(lowestPenality * distance);
    weights[secondNode].push_back(lowestPenality * distance);
  }
  return bellmanFordShortestPath(graph, weights, start, end);
}

std::vector<std::vector<size_t>>
RegionDivider::computePaths(std::vector<NodeId> corners) {
  assert(corners.size() == nbRegion);
  auto sumBools = [](std::array<bool, DIMENSION> &v) {
    return v[0] + v[1] + v[2];
  };
  std::vector<std::vector<size_t>> paths;
  for (RegionId startRegionId = 0; startRegionId < nbRegion; startRegionId++) {
    std::vector<std::array<bool, DIMENSION>> sharedLoops(nbRegion,
                                                         {false, false, false});
    for (Axis d = 0; d < DIMENSION; d++) {
      ZoneId zoneId = loopZones[d].getZone(startRegionId);
      for (RegionId other : loopZones[d][zoneId])
        sharedLoops[other][d] = true;
    }
    for (RegionId endRegionId = 0; endRegionId < nbRegion; endRegionId++) {
      if (endRegionId == startRegionId)
        continue;
      size_t score = sumBools(sharedLoops[endRegionId]);
      assert(score < 3);
      if (score != 2 || endRegionId > startRegionId)
        continue;
      Axis a = 0;
      while (sharedLoops[endRegionId][a])
        a++;
      std::vector<size_t> path =
          linkCorners(corners[startRegionId], corners[endRegionId], a);
      paths.push_back(path);
    }
  }

  std::cout << "GOT " << paths.size() << " PATHS" << std::endl;
  return paths;
}

void RegionDivider::saveMesh(std::string filename) {
  std::vector<NodeId> corners = computeCorners();

  std::vector<size_t> nodePathId(mesh.vertices.size(), 0);
  std::vector<std::vector<size_t>> paths = computePaths(corners);
  savePaths("primary_paths.vtu", paths, mesh.vertices);

  std::vector<Vertex> centerEdges(mesh.edges.size());
  for (size_t i = 0; i < mesh.edges.size(); i++)
    centerEdges[i] = mesh.edges[i].computeCenter(mesh);

  std::vector<bool> isCorner(mesh.vertices.size(), false);

  std::vector<ZoneId> xZone(mesh.vertices.size(), false);
  std::vector<int> xSign(mesh.vertices.size(), NAN);
  std::vector<ZoneId> yZone(mesh.vertices.size(), false);
  std::vector<int> ySign(mesh.vertices.size(), NAN);
  std::vector<ZoneId> zZone(mesh.vertices.size(), false);
  std::vector<int> zSign(mesh.vertices.size(), NAN);
  for (NodeId nodeId = 0; nodeId < mesh.nbBoundaryNode; nodeId++) {
    RegionId regionId = regions[nodeId];
    xZone[nodeId] = loopZones[0].getZone(regionId);
    yZone[nodeId] = loopZones[1].getZone(regionId);
    zZone[nodeId] = loopZones[2].getZone(regionId);

    xSign[nodeId] = getSignedAxisValue(regionToNormals[regionId][0]);
    ySign[nodeId] = getSignedAxisValue(regionToNormals[regionId][1]);
    zSign[nodeId] = getSignedAxisValue(regionToNormals[regionId][2]);
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

  saver.addField<ElementType::Node>("xSign", xSign);
  saver.addField<ElementType::Node>("ySign", ySign);
  saver.addField<ElementType::Node>("zSign", zSign);

  saver.addField<ElementType::Node>("pathId", nodePathId);
  saver.saveMesh();
  saver.saveEdges();
}
