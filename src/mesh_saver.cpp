#include "mesh_saver.hpp"
#include <vtkIdTypeArray.h>

std::vector<Edge> buildEdges(const std::vector<Tetrahedron> &tets) {
  std::vector<Edge> edges;

  for (const auto &t : tets) {
    for (size_t i = 0; i < NB_NODE_PER_TET; ++i) {
      for (size_t j = i + 1; j < NB_NODE_PER_TET; ++j) {
        Edge e;
        e.nodes[0] = std::min(t[i], t[j]);
        e.nodes[1] = std::max(t[i], t[j]);
        edges.push_back(e);
      }
    }
  }
  return edges;
}

MeshSaver::MeshSaver(const std::vector<Vertex> &vertices,
                     const std::vector<Edge> &edges,
                     const std::vector<Tetrahedron> &tets,
                     const std::string &path)
    : path(path), nbNodes(vertices.size()), nbCells(tets.size()) {
  grid = vtkSmartPointer<vtkUnstructuredGrid>::New();

  auto points = vtkSmartPointer<vtkPoints>::New();

  for (const auto &v : vertices) {
    points->InsertNextPoint(v.x(), v.y(), v.z());
  }

  grid->SetPoints(points);

  for (const auto &t : tets) {
    vtkIdType ids[NB_NODE_PER_TET];

    for (size_t i = 0; i < NB_NODE_PER_TET; ++i) {
      if (t[i] >= vertices.size()) {
        throw std::out_of_range("Invalid tetrahedron vertex index.");
      }

      ids[i] = static_cast<vtkIdType>(t[i]);
    }

    grid->InsertNextCell(VTK_TETRA, NB_NODE_PER_TET, ids);
  }

  nbEdges = edges.size();

  edgeGrid = vtkSmartPointer<vtkPolyData>::New();
  edgeGrid->SetPoints(points);

  auto lines = vtkSmartPointer<vtkCellArray>::New();

  for (const auto &edge : edges) {
    vtkIdType ids[2] = {static_cast<vtkIdType>(edge[0]),
                        static_cast<vtkIdType>(edge[1])};

    lines->InsertNextCell(2, ids);
  }

  edgeGrid->SetLines(lines);
}

void MeshSaver::saveMesh() const {
  auto writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();

  writer->SetFileName(path.c_str());
  writer->SetInputData(grid);
  writer->SetDataModeToBinary();

  if (writer->Write() != 1) {
    throw std::runtime_error("Failed to write VTK mesh to " + path);
  }
}

void MeshSaver::saveEdges() const {
  assert(nbEdges > 0);
  std::string edgePath = path;

  const auto pos = edgePath.rfind(".vtu");

  if (pos != std::string::npos) {
    edgePath.replace(pos, 4, "_edges.vtp");
  } else {
    edgePath += "_edges.vtp";
  }

  auto edgeWriter = vtkSmartPointer<vtkXMLPolyDataWriter>::New();

  edgeWriter->SetFileName(edgePath.c_str());
  edgeWriter->SetInputData(edgeGrid);
  edgeWriter->SetDataModeToBinary();

  if (edgeWriter->Write() != 1) {
    throw std::runtime_error("Failed to write VTK edges to " + edgePath);
  }
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
  axes->SetName("Axis");
  axes->SetNumberOfComponents(1);

  // Insert paths as polyline cells.
  for (const auto &[nodes, axis] : paths) {
    if (nodes.size() < 3)
      throw std::runtime_error("A closed path must contain at least 3 nodes.");

    if (nodes.front() != nodes.back())
      throw std::runtime_error("Path is not closed.");

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
