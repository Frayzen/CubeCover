#include "span_tree.hpp"
#include "constants.hpp"
#include <stack>

void convertToSpanTree(TetrahedralMesh &mesh) {

  std::vector<bool> tetsPushed(mesh.tets.size());
  std::fill(tetsPushed.begin(), tetsPushed.end(), false);
  std::stack<size_t> tetsToProcess;
  tetsToProcess.push(0);
  tetsPushed[0] = true;

  while (!tetsToProcess.empty()) {
    size_t curTet = tetsToProcess.top();
    tetsToProcess.pop();

    for (auto face : mesh.tetToFaces[curTet]) {
      Face &f = mesh.faces[face];
      size_t parentTet = f.adjacentCells[0];
      if (parentTet != curTet) {
        if (!tetsPushed[parentTet])
          std::swap(f.adjacentCells[0], f.adjacentCells[1]);
      }
      size_t childTet = f.adjacentCells[1];
      if (childTet != curTet && childTet != UNDEFINED_ELEMENT &&
          !tetsPushed[childTet]) {
        tetsToProcess.push(childTet);
        tetsPushed[childTet] = true;
      }
      tetsPushed[curTet] = true;
    }
  }
}
