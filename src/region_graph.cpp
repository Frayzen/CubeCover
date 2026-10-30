#include "region_graph.hpp"
#include "mesh_saver.hpp"
#include <queue>
#include <string>
RegionGraph::RegionGraph(EdgeGraph &g)
    : g(g), mesh(g.mesh), nbRegion(1), regions() {
  for (size_t i = 0; i < g.mesh.nbBoundaryNode; i++) {
    regions.set(i, 0);
  }
}

void RegionGraph::fillConnectedFaces(NodeId startingEdgeId,
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
  for (auto [adjEdgeId, _] : g.adjacency[startingEdgeId]) {
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
    for (auto [adjEdgeId, _] : g.adjacency[edgeId]) {
      if (processed[adjEdgeId])
        continue;
      if (mesh.edges[adjEdgeId].containsNode(followUpNode)) {
        toProcess.push(adjEdgeId);
        processed[adjEdgeId] = true;
      }
    }
  }
}

void RegionGraph::insertLoop(std::vector<size_t> edgeIds, size_t axis) {
  paths.push_back(std::make_pair(edgeIds, axis));

  for (size_t i = 0; i < edgeIds.size() - 1; i++) {
    size_t edgeId = edgeIds[i];

    Edge &e = mesh.edges[edgeId];
    bool sameRegion = regions[e[0]] == regions[e[1]];
    if (!sameRegion)
      continue;
    fillConnectedFaces(edgeId, edgeIds, axis);
    std::cout << "GOT " << nbRegion << " REGIONS " << std::endl;
  }
}

void RegionGraph::saveMesh(std::string filename) {

  std::vector<int> colors(mesh.vertices.size(), -1);
  for (size_t nodeId = 0; nodeId < mesh.nbBoundaryNode; nodeId++)
    colors[nodeId] = static_cast<size_t>(getFaceRegion(nodeId));

  MeshSaver saver(mesh.vertices, mesh.edges, mesh.tets, "regions.vtu");
  saver.addField<ElementType::Node>("region", colors);
  saver.saveMesh();
}
