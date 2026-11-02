#include "graph.hpp"
#include "region_divider.hpp"
#include <queue>

using AxisFlag = std::array<bool, DIMENSION>;

std::vector<PerDimNormals>
RegionDivider::colorizeMesh(std::vector<Path> &paths) {
  Graph nodeGraph(mesh.nbBoundaryNode);
  for (EdgeId edgeId = 0; edgeId < mesh.nbBoundaryEdge; edgeId++) {
    Edge &e = mesh.edges[edgeId];
    nodeGraph.adjacency[e[0]].push_back(e[1]);
    nodeGraph.adjacency[e[1]].push_back(e[0]);
  }
  std::vector<AxisFlag> nodeToAxisFlag(mesh.nbBoundaryNode, {false});
  std::vector<bool> processed(mesh.nbBoundaryNode, false);
  for (auto [nodes, axis] : paths) {
    for (NodeId nodeId : nodes) {
      processed[nodeId] = true;
      for (Axis d = 0; d < DIMENSION; d++) {
        if (d != axis)
          nodeToAxisFlag[nodeId][axis] = true;
      }
    }
  }
  auto isCorner = [&](NodeId nodeId) {
    return processed[nodeId] && nodeToAxisFlag[nodeId][0] &&
           nodeToAxisFlag[nodeId][1] && nodeToAxisFlag[nodeId][2];
  };
  std::queue<NodeId> toProcess;
  for (NodeId initNodeId = 0; initNodeId < mesh.nbBoundaryNode; initNodeId++) {
    if (processed[initNodeId])
      continue;
    toProcess.push(initNodeId);
    processed[initNodeId] = true;
    AxisFlag batchAxisFlag = {true, true, true};
    std::vector<NodeId> batch;
    while (!toProcess.empty()) {
      NodeId nodeId = toProcess.front();
      batch.push_back(nodeId);
      toProcess.pop();
      for (NodeId adjId : nodeGraph.adjacency[nodeId]) {
        if (isCorner(adjId))
          continue;
        for (Axis d = 0; d < DIMENSION; d++)
          batchAxisFlag[d] &= !nodeToAxisFlag[adjId][d];
        if (processed[adjId])
          continue;
        toProcess.push(adjId);
        processed[adjId] = true;
      }
    }
    for (NodeId nodeId : batch)
      nodeToAxisFlag[nodeId] = batchAxisFlag;
  }
  std::vector<PerDimNormals> res(mesh.nbBoundaryNode, {ZERO, ZERO, ZERO});
  for (NodeId nodeId = 0; nodeId < mesh.nbBoundaryNode; nodeId++) {
    for (Axis d = 0; d < DIMENSION; d++) {
      if (!nodeToAxisFlag[nodeId][d])
        continue;
      res[nodeId][d] = regionToNormals[regions[nodeId]][d];
    }
  }
  for (auto [nodes, axis] : paths) {
    for (NodeId nodeId : nodes) {
      res[nodeId] = {ZERO, ZERO, ZERO};
    }
  }
  for (auto [nodes, axis] : paths) {
    for (NodeId nodeId : nodes) {
      for (Axis d = 0; d < DIMENSION; d++)
        if (d != axis)
          res[nodeId][d] = regionToNormals[regions[nodeId]][d];
    }
  }

  return res;
}
