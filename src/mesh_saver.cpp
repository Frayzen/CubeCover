#include "mesh_saver.hpp"
#include "mesh_loader.hpp"
#include "vertex.hpp"
#include <stdexcept>

#include <vtkIntArray.h>
#include <vtkPointData.h>

#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkCellType.h>
#include <vtkDoubleArray.h>
#include <vtkIntArray.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkSmartPointer.h>
#include <vtkTetra.h>
#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridWriter.h>
#include <vtkXMLUnstructuredGridWriter.h>

#include <Eigen/Dense>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

void saveMeshVtk(const std::vector<Vertex> &vertices,
                 const std::vector<Tetrahedron> &tets,
                 const std::vector<Eigen::Matrix3d> &frames,
                 const std::string &path) {
  if (frames.size() != tets.size()) {
    throw std::invalid_argument(
        "The number of frames must match the number of tetrahedra.");
  }

  auto points = vtkSmartPointer<vtkPoints>::New();
  auto grid = vtkSmartPointer<vtkUnstructuredGrid>::New();

  // Add unique mesh vertices.
  for (const auto &v : vertices) {
    points->InsertNextPoint(v.x(), v.y(), v.z());
  }

  grid->SetPoints(points);

  // Add tetrahedra using mesh connectivity.
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

  // Create frame arrays.
  std::array<vtkSmartPointer<vtkDoubleArray>, 3> frameArrays;

  for (size_t j = 0; j < 3; ++j) {
    frameArrays[j] = vtkSmartPointer<vtkDoubleArray>::New();
    frameArrays[j]->SetName(("frame_col" + std::to_string(j)).c_str());
    frameArrays[j]->SetNumberOfComponents(3);
    frameArrays[j]->SetNumberOfTuples(static_cast<vtkIdType>(frames.size()));
  }

  for (vtkIdType i = 0; i < static_cast<vtkIdType>(frames.size()); ++i) {
    const Eigen::Matrix3d &F = frames[i];

    for (size_t j = 0; j < 3; ++j) {
      const double v[3] = {F(0, j), F(1, j), F(2, j)};

      frameArrays[j]->SetTuple(i, v);
    }
  }

  // Store frames as cell data.
  for (auto &array : frameArrays) {
    grid->GetCellData()->AddArray(array);
  }

  auto writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();

  writer->SetFileName(path.c_str());
  writer->SetInputData(grid);
  writer->SetDataModeToBinary();

  if (writer->Write() != 1) {
    throw std::runtime_error("Failed to write VTK mesh to " + path);
  }
}
