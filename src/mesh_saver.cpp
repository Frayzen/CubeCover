#include "mesh_saver.hpp"

#include <cassert>
#include <vector>

#include <vtkIdTypeArray.h>

MeshSaver::MeshSaver(const TetrahedralMesh &mesh, const std::string &path)
    : mesh(mesh), path(path) {

  buildVolumeGrid();
  buildEdgeGrid();
  buildFaceGrid();
}

void MeshSaver::buildVolumeGrid() {

  grid = vtkSmartPointer<vtkUnstructuredGrid>::New();

  auto points = vtkSmartPointer<vtkPoints>::New();

  for (const auto &v : mesh.vertices) {
    points->InsertNextPoint(v.x(), v.y(), v.z());
  }

  grid->SetPoints(points);

  for (const auto &tet : mesh.tets) {

    vtkIdType ids[NB_NODE_PER_TET];

    for (size_t i = 0; i < NB_NODE_PER_TET; ++i) {

      if (tet[i] >= mesh.vertices.size()) {
        throw std::out_of_range("Invalid tetrahedron vertex index.");
      }

      ids[i] = static_cast<vtkIdType>(tet[i]);
    }

    grid->InsertNextCell(VTK_TETRA, NB_NODE_PER_TET, ids);
  }
}

void MeshSaver::buildEdgeGrid() {

  edgeGrid = vtkSmartPointer<vtkPolyData>::New();
  edgeGrid->SetPoints(grid->GetPoints());

  auto lines = vtkSmartPointer<vtkCellArray>::New();

  for (size_t i = 0; i < mesh.nbBoundaryEdge; ++i) {

    const auto &edge = mesh.edges[i];

    if (edge[0] >= mesh.vertices.size() || edge[1] >= mesh.vertices.size()) {
      throw std::out_of_range("Invalid edge vertex index.");
    }

    vtkIdType ids[2] = {static_cast<vtkIdType>(edge[0]),
                        static_cast<vtkIdType>(edge[1])};

    lines->InsertNextCell(2, ids);
  }

  edgeGrid->SetLines(lines);
}

void MeshSaver::buildFaceGrid() {

  faceGrid = vtkSmartPointer<vtkPolyData>::New();
  faceGrid->SetPoints(grid->GetPoints());

  auto triangles = vtkSmartPointer<vtkCellArray>::New();

  for (size_t i = 0; i < mesh.nbBoundaryFace; ++i) {

    const auto &face = mesh.faces[i];

    for (size_t j = 0; j < 3; ++j) {
      if (face[j] >= mesh.vertices.size()) {
        throw std::out_of_range("Invalid face vertex index.");
      }
    }

    vtkIdType ids[3] = {static_cast<vtkIdType>(face[0]),
                        static_cast<vtkIdType>(face[1]),
                        static_cast<vtkIdType>(face[2])};

    triangles->InsertNextCell(3, ids);
  }

  faceGrid->SetPolys(triangles);
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

  assert(mesh.nbBoundaryEdge > 0);

  std::string edgePath = path;

  const auto pos = edgePath.rfind(".vtu");

  if (pos != std::string::npos) {
    edgePath.replace(pos, 4, "_edges.vtp");
  } else {
    edgePath += "_edges.vtp";
  }

  auto writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();

  writer->SetFileName(edgePath.c_str());
  writer->SetInputData(edgeGrid);
  writer->SetDataModeToBinary();

  if (writer->Write() != 1) {
    throw std::runtime_error("Failed to write VTK edges to " + edgePath);
  }
}

void MeshSaver::saveFaces() const {

  assert(mesh.nbBoundaryFace > 0);

  std::string facePath = path;

  const auto pos = facePath.rfind(".vtu");

  if (pos != std::string::npos) {
    facePath.replace(pos, 4, "_faces.vtp");
  } else {
    facePath += "_faces.vtp";
  }

  auto writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();

  writer->SetFileName(facePath.c_str());
  writer->SetInputData(faceGrid);
  writer->SetDataModeToBinary();

  if (writer->Write() != 1) {
    throw std::runtime_error("Failed to write VTK faces to " + facePath);
  }
}
