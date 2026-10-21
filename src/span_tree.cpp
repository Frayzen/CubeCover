
#include "span_tree.hpp"
#include "constants.hpp"
#include "mesh_loader.hpp"

#include <stack>
#include <stdexcept>
#include <vector>

std::vector<size_t> convertToSpanTree(TetrahedralMesh &mesh) {
  const size_t n = mesh.tets.size();
  std::vector<bool> visited(n, false);
  std::vector<size_t> spanTree(mesh.amountInternalFace);
  size_t curSpanTreeFace = 0;
  std::stack<size_t> pending;
  const size_t root = 373;
  pending.push(root);
  visited[root] = true;

  while (!pending.empty()) {
    const size_t curTet = pending.top();
    pending.pop();

    for (size_t faceId : mesh.tetToFaces[curTet]) {
      Face &face = mesh.faces[faceId];

      size_t a = face.adjacentCells[0];
      size_t b = face.adjacentCells[1];
      if (b == UNDEFINED_ELEMENT) // Boundary face
        continue;

      size_t neighbor = a == curTet ? b : a;
      if (neighbor == curTet)
        throw std::runtime_error("Inconsistent tetToFaces adjacency");
      if (visited[neighbor])
        continue;

      spanTree[curSpanTreeFace++] = faceId;
      face.adjacentCells[0] = curTet;
      face.adjacentCells[1] = neighbor;

      visited[neighbor] = true;
      pending.push(neighbor);
    }
  }
  spanTree.resize(curSpanTreeFace);
  return spanTree;
}
