#include "path_saver.hpp"
#include <vtkCellData.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLUnstructuredGridWriter.h>

void savePaths(std::string filename, std::vector<std::vector<size_t>> paths,
               std::vector<Vertex> vertices) {
  std::vector<std::pair<std::vector<NodeId>, Axis>> pathUpdated;
  for (size_t i = 0; i < paths.size(); i++) {
    pathUpdated.push_back(std::make_pair(paths[i], i));
  }
  return savePaths(filename, pathUpdated, vertices);
}

void savePaths(std::string filename,
               std::vector<std::pair<std::vector<NodeId>, Axis>> paths,
               std::vector<Vertex> vertices) {
  auto points = vtkSmartPointer<vtkPoints>::New();
  points->SetNumberOfPoints(vertices.size());

  // Insert vertices.
  for (size_t i = 0; i < vertices.size(); ++i) {
    points->SetPoint(i, vertices[i][0], vertices[i][1], vertices[i][2]);
  }

  auto grid = vtkSmartPointer<vtkUnstructuredGrid>::New();
  grid->SetPoints(points);

  // Axis associated with each path (cell data).
  auto axes = vtkSmartPointer<vtkIdTypeArray>::New();
  axes->SetName("Id");
  axes->SetNumberOfComponents(1);

  // Insert paths as polyline cells.
  for (const auto &[nodes, axis] : paths) {
    if (nodes.size() < 3)
      throw std::runtime_error("A path must contain at least 3 nodes.");

    std::vector<vtkIdType> ids;
    ids.reserve(nodes.size());

    for (NodeId node : nodes) {
      if (node >= vertices.size())
        throw std::out_of_range("Path references an invalid vertex.");

      ids.push_back(static_cast<vtkIdType>(node));
    }

    grid->InsertNextCell(VTK_POLY_LINE, static_cast<vtkIdType>(ids.size()),
                         ids.data());

    axes->InsertNextValue(static_cast<vtkIdType>(axis));
  }

  grid->GetCellData()->AddArray(axes);

  // Write VTU file.
  auto writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
  writer->SetFileName(filename.c_str());
  writer->SetInputData(grid);

  if (writer->Write() != 1)
    throw std::runtime_error("Failed to write VTU file: " + filename);
}
