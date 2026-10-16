#include "mesh_saver.hpp"
#include <stdexcept>

#include <vtkCellArray.h>
#include <vtkCellType.h>
#include <vtkPoints.h>
#include <vtkSmartPointer.h>
#include <vtkTetra.h>
#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridWriter.h>

#include <Eigen/Dense>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

void saveMesh(const std::vector<Eigen::Vector3d> &vertices,
              const std::string &filename) {
  if (vertices.size() % 4 != 0) {
    throw std::invalid_argument(
        "The number of vertices must be a multiple of 4.");
  }

  std::ofstream out(filename);
  if (!out) {
    throw std::runtime_error("Cannot open file: " + filename);
  }

  const std::size_t numVertices = vertices.size();
  const std::size_t numTets = numVertices / 4;

  // Medit ASCII mesh format
  out << "MeshVersionFormatted 2\n";
  out << "Dimension 3\n\n";

  out << "Vertices\n" << numVertices << '\n';

  for (const auto &v : vertices) {
    out << v.x() << ' ' << v.y() << ' ' << v.z() << " 0\n";
  }

  out << "\nTetrahedra\n" << numTets << '\n';

  // Medit uses 1-based vertex indices.
  for (std::size_t i = 0; i < numTets; ++i) {
    const std::size_t base = 4 * i + 1;
    out << base << ' ' << base + 1 << ' ' << base + 2 << ' ' << base + 3
        << " 0\n";
  }

  out << "\nEnd\n";

  if (!out) {
    throw std::runtime_error("Failed to write mesh: " + filename);
  }
}
