#include "surface_mesh.hpp"
#include "voxel_mesh.hpp"

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <functional>
#include <queue>
#include <stdexcept>
#include <utility>
#include <vector>

SurfaceMesh::SurfaceMesh(const VoxelMesh &mesh)
    : m(mesh), regionToQuads(m.regions.size()) {}

// ============================================================
// Find region corner
// ============================================================

NodeId SurfaceMesh::findCorner(RegionId regionId) const {

  const Region &r = m.regions[regionId];

  Vertex target = Vertex::Zero();

  int numNormals = 0;

  for (Axis d = 0; d < DIMENSION; ++d) {

    assert(r.signedNormals[d] != UNDEFINED);

    if (r.signedNormals[d] == ZERO) {

      target[d] = (r.startSlices[d] + r.endSlices[d]) / 2.0;

      continue;
    }

    ++numNormals;

    target[d] =
        r.signedNormals[d] == POSITIVE ? r.endSlices[d] : r.startSlices[d];
  }
  target *= m.voxelSize;

  assert(numNormals > 0);

  double minDist = std::numeric_limits<double>::max();
  NodeId minVertex = NO_NODE;

  for (QuadId quadId : regionToQuads[regionId]) {
    for (EdgeId edgeId : quadToEdges[quadId]) {
      for (NodeId vertexId : edges[edgeId]) {
        const Vertex &v = vertices[vertexId];
        const double dist = (v - target).norm();

        if (dist < minDist) {
          minDist = dist;
          minVertex = vertexId;
        }
      }
    }
  }

  assert(minVertex != NO_NODE);

  return minVertex;
}

// ============================================================
// Find all region corners
// ============================================================

std::vector<NodeId> SurfaceMesh::findCorners() const {

  std::vector<NodeId> res;

  res.reserve(m.regions.size());

  for (RegionId rId = 0; rId < m.regions.size(); ++rId)
    res.push_back(findCorner(rId));

  return res;
}

std::vector<NodeId> SurfaceMesh::getNeighbours(NodeId vertexId) const {

  std::vector<NodeId> res;

  for (EdgeId eId = 0; eId < edges.size(); ++eId) {
    const auto &edge = edges[eId];

    if (edge[0] == vertexId)
      res.push_back(edge[1]);
    else if (edge[1] == vertexId)
      res.push_back(edge[0]);
  }

  return res;
}

// ============================================================
// Dijkstra shortest path
// ============================================================

std::vector<EdgeId> SurfaceMesh::findShortestPath(
    NodeId start, NodeId goal,
    std::array<double, DIMENSION> axisPenality) const {

  using QueueItem = std::pair<double, NodeId>;

  const double INF = std::numeric_limits<double>::max();

  std::vector<double> dist(vertices.size(), INF);

  std::vector<NodeId> previous(vertices.size(), NO_NODE);

  std::vector<EdgeId> previousEdge(vertices.size(), NO_EDGE);

  std::priority_queue<QueueItem, std::vector<QueueItem>,
                      std::greater<QueueItem>>
      queue;

  dist[start] = 0.0;

  queue.emplace(0.0, start);

  while (!queue.empty()) {

    const auto [currentDist, current] = queue.top();

    queue.pop();

    if (currentDist > dist[current])
      continue;

    if (current == goal)
      break;

    for (EdgeId eId : vertexToEdges[current]) {

      const Edge &edge = edges[eId];

      const NodeId neighbor = edge[0] == current ? edge[1] : edge[0];

      Vertex dir = vertices[current] - vertices[neighbor];
      double weight = dir.norm();

      std::array<bool, DIMENSION> foundAxis = {false, false, false};
      size_t amountDifferentAxis = 0;
      Axis dirAxis = NO_AXIS;
      for (QuadId qId : edgeToQuads[eId]) {
        Axis a = quadToNormal[qId];
        if (!foundAxis[a])
          amountDifferentAxis++;
        if (dir[a] > m.voxelSize / 2)
          dirAxis = a;
        foundAxis[a] = true;
      }

      // Done in order to prefer sharp edges rather than planar surfaces
      if (amountDifferentAxis == 1)
        weight *= (1 + weight / 100.0);
      weight *= axisPenality[dirAxis];

      const double candidate = currentDist + weight;

      if (candidate < dist[neighbor]) {

        dist[neighbor] = candidate;

        previous[neighbor] = current;
        previousEdge[neighbor] = eId;

        queue.emplace(candidate, neighbor);
      }
    }
  }

  if (start != goal && previous[goal] == NO_NODE)
    return {};

  std::vector<EdgeId> path;

  for (NodeId current = goal; current != start; current = previous[current]) {

    if (previousEdge[current] == NO_EDGE)
      return {};

    path.push_back(previousEdge[current]);
  }

  std::reverse(path.begin(), path.end());

  return path;
}

// ============================================================
// Find paths between region corners
// ============================================================

std::vector<SurfacePath> SurfaceMesh::findPaths() {

  regionCorners = findCorners();

  paths.clear();

  for (RegionId rId = 0; rId < m.regions.size(); ++rId) {

    const Region &r = m.regions[rId];

    for (Axis a = 0; a < DIMENSION; ++a) {

      for (const RegionInterface &interface : r.children[a]) {

        const RegionId adjId = interface.regionId;
        const Region &adj = m.regions[adjId];

        SurfacePath path;
        path.regions[0] = rId;
        path.regions[1] = adjId;

        path.normal = interface.normal;

        std::array<double, DIMENSION> axisPenality;
        for (Axis a = 0; a < DIMENSION; a++) {
          if (a == interface.normal)
            axisPenality[a] = 0;
          else if (adj.signedNormals[a] != r.signedNormals[a])
            axisPenality[a] = 1;
          else
            axisPenality[a] = 2;
        }

        path.edges = findShortestPath(regionCorners[rId], regionCorners[adjId],
                                      axisPenality);

        const PathId pathId = paths.size();

        for (EdgeId eId : path.edges)
          edgeToPaths[eId].push_back(pathId);

        paths.push_back(std::move(path));
      }
    }
  }

  return paths;
}

std::vector<Axis> SurfaceMesh::colorize() const {

  if (quadToEdges.size() != quads.size())
    throw std::runtime_error("SurfaceMesh: quadToEdges size mismatch.");

  if (edgeToQuads.size() != edges.size())
    throw std::runtime_error("SurfaceMesh: edgeToQuads size mismatch.");

  if (edgeToPaths.size() != edges.size())
    throw std::runtime_error("SurfaceMesh: edgeToPaths size mismatch.");

  std::vector<bool> visited(quads.size(), false);

  std::vector<Axis> colors(quads.size(), NO_AXIS);

  for (QuadId seed = 0; seed < quads.size(); ++seed) {

    if (visited[seed])
      continue;

    // Axes of paths bordering this component.
    std::array<size_t, DIMENSION> encountered{};

    std::vector<QuadId> component;
    std::queue<QuadId> queue;

    visited[seed] = true;
    queue.push(seed);
    RegionId seedRegionId = quadToRegion[seed];

    // ----------------------------------------------------
    // Flood fill
    // ----------------------------------------------------

    while (!queue.empty()) {

      const QuadId qId = queue.front();
      queue.pop();

      component.push_back(qId);

      // Inspect the four edges of the current quad.
      for (EdgeId eId : quadToEdges[qId]) {

        if (eId >= edges.size())
          throw std::runtime_error("SurfaceMesh: invalid quad edge index.");

        // ------------------------------------------------
        // Path edge: record its normal and do not cross.
        // ------------------------------------------------

        if (!edgeToPaths[eId].empty()) {

          for (PathId pId : edgeToPaths[eId]) {
            if (pId >= paths.size())
              throw std::runtime_error("SurfaceMesh: invalid path index.");

            if (paths[pId].regions[0] != seedRegionId &&
                paths[pId].regions[1] != seedRegionId)
              continue;

            const Axis a = paths[pId].normal;
            const std::size_t axis = static_cast<std::size_t>(a);

            if (axis >= DIMENSION)
              throw std::runtime_error("SurfaceMesh: invalid path normal.");

            encountered[axis]++;
          }

          continue;
        }

        // ------------------------------------------------
        // Non-path edge: visit adjacent quads.
        // ------------------------------------------------

        for (QuadId neighbor : edgeToQuads[eId]) {

          if (neighbor >= quads.size())
            throw std::runtime_error(
                "SurfaceMesh: invalid adjacent quad index.");

          if (visited[neighbor])
            continue;

          visited[neighbor] = true;
          queue.push(neighbor);
        }
      }
    }

    // ----------------------------------------------------
    // Find the missing axis
    // ----------------------------------------------------

    // TODO improve this heuristic
    Axis missingAxis = NO_AXIS;
    size_t minimumCounts = std::numeric_limits<size_t>::max();

    for (Axis a = 0; a < DIMENSION; ++a) {

      if (encountered[a] >= minimumCounts)
        continue;
      minimumCounts = encountered[a];
      missingAxis = a;
    }

    // ----------------------------------------------------
    // Assign the missing axis to every quad in the component
    // ----------------------------------------------------

    for (QuadId qId : component)
      colors[qId] = missingAxis;
  }

  return colors;
}
