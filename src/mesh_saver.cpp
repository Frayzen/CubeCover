#include "mesh_saver.hpp"
#include "constants.hpp"
#include <stdexcept>

#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkCellType.h>
#include <vtkDoubleArray.h>
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

void saveMeshGmsh(const std::vector<Eigen::Vector3d> &vertices,
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

#include <stdexcept>
#include <string>
#include <vector>

#include <Eigen/Dense>

#include <vtkCellData.h>
#include <vtkDoubleArray.h>
#include <vtkPoints.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridWriter.h>

void saveMeshVtk(const std::vector<Eigen::Vector3d> &vertices,
                 const std::vector<Eigen::Matrix3d> &frames,
                 const std::string &path) {
  if (vertices.size() % 4 != 0) {
    throw std::invalid_argument(
        "The number of vertices must be a multiple of 4.");
  }

  const std::size_t numTetrahedra = vertices.size() / 4;

  if (frames.size() != numTetrahedra) {
    throw std::invalid_argument(
        "The number of frames must match the number of tetrahedra.");
  }

  auto points = vtkSmartPointer<vtkPoints>::New();
  auto grid = vtkSmartPointer<vtkUnstructuredGrid>::New();

  // Add vertices.
  for (const auto &v : vertices) {
    points->InsertNextPoint(v.x(), v.y(), v.z());
  }

  grid->SetPoints(points);

  // Add tetrahedra.
  for (vtkIdType i = 0; i < static_cast<vtkIdType>(vertices.size()); i += 4) {
    vtkIdType ids[4] = {i, i + 1, i + 2, i + 3};
    grid->InsertNextCell(VTK_TETRA, 4, ids);
  }

  // Create three vector arrays, one for each frame column.
  auto col0 = vtkSmartPointer<vtkDoubleArray>::New();
  auto col1 = vtkSmartPointer<vtkDoubleArray>::New();
  auto col2 = vtkSmartPointer<vtkDoubleArray>::New();

  col0->SetName("frame_col0");
  col1->SetName("frame_col1");
  col2->SetName("frame_col2");

  col0->SetNumberOfComponents(3);
  col1->SetNumberOfComponents(3);
  col2->SetNumberOfComponents(3);

  col0->SetNumberOfTuples(static_cast<vtkIdType>(frames.size()));
  col1->SetNumberOfTuples(static_cast<vtkIdType>(frames.size()));
  col2->SetNumberOfTuples(static_cast<vtkIdType>(frames.size()));

  for (vtkIdType i = 0; i < static_cast<vtkIdType>(frames.size()); ++i) {
    const Eigen::Matrix3d &F = frames[i];

    const double v0[3] = {F(0, 0), F(1, 0), F(2, 0)};
    const double v1[3] = {F(0, 1), F(1, 1), F(2, 1)};
    const double v2[3] = {F(0, 2), F(1, 2), F(2, 2)};

    col0->SetTuple(i, v0);
    col1->SetTuple(i, v1);
    col2->SetTuple(i, v2);
  }

  // Store vectors as cell data.
  grid->GetCellData()->AddArray(col0);
  grid->GetCellData()->AddArray(col1);
  grid->GetCellData()->AddArray(col2);

  // Write VTK mesh.
  auto writer = vtkSmartPointer<vtkUnstructuredGridWriter>::New();
  writer->SetFileName(path.c_str());
  writer->SetInputData(grid);
  writer->SetFileTypeToBinary();

  if (writer->Write() != 1) {
    throw std::runtime_error("Failed to write VTK mesh to " + path);
  }
}
