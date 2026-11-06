#ifndef VOXEL_MESH_SAVER_HPP
#define VOXEL_MESH_SAVER_HPP

#include "voxel_mesh.hpp"

#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

#include <Eigen/Dense>

#include <vtkCellData.h>
#include <vtkDataSetAttributes.h>
#include <vtkDoubleArray.h>
#include <vtkPoints.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLUnstructuredGridWriter.h>

class VoxelMeshSaver {
private:
  const VoxelMesh &mesh;

  vtkSmartPointer<vtkUnstructuredGrid> grid;

  std::string path;

  template <typename T> static constexpr int nbComponents() {

    if constexpr (std::is_arithmetic_v<T>) {

      return 1;

    } else if constexpr (requires {
                           T::RowsAtCompileTime;
                           T::ColsAtCompileTime;
                         }) {

      static_assert(T::RowsAtCompileTime != Eigen::Dynamic &&
                        T::ColsAtCompileTime != Eigen::Dynamic,
                    "Dynamic-size Eigen types are not supported.");

      return T::RowsAtCompileTime * T::ColsAtCompileTime;

    } else {

      static_assert(sizeof(T) == 0, "Unsupported field type");
    }
  }

  template <typename T> static double component(const T &value, size_t i) {

    if constexpr (std::is_arithmetic_v<T>) {

      return static_cast<double>(value);

    } else {

      constexpr size_t rows = T::RowsAtCompileTime;

      return static_cast<double>(value(i % rows, i / rows));
    }
  }

  template <typename T>
  void addFieldToArray(vtkDataSetAttributes *attributes,
                       const std::string &name, const std::vector<T> &values) {

    auto array = vtkSmartPointer<vtkDoubleArray>::New();

    constexpr int nComp = nbComponents<T>();

    array->SetName(name.c_str());
    array->SetNumberOfComponents(nComp);
    array->SetNumberOfTuples(static_cast<vtkIdType>(values.size()));

    for (size_t i = 0; i < values.size(); ++i) {

      double tuple[nComp];

      for (int j = 0; j < nComp; ++j) {
        tuple[j] = component(values[i], j);
      }

      array->SetTuple(static_cast<vtkIdType>(i), tuple);
    }

    attributes->AddArray(array);
  }

  void buildGrid() {

    grid = vtkSmartPointer<vtkUnstructuredGrid>::New();

    auto points = vtkSmartPointer<vtkPoints>::New();

    const double s = mesh.voxelSize;

    for (const auto &voxel : mesh.voxels) {

      const double x = voxel[0] * s;
      const double y = voxel[1] * s;
      const double z = voxel[2] * s;

      // VTK_HEXAHEDRON vertex ordering

      points->InsertNextPoint(x, y, z);
      points->InsertNextPoint(x + s, y, z);
      points->InsertNextPoint(x + s, y + s, z);
      points->InsertNextPoint(x, y + s, z);

      points->InsertNextPoint(x, y, z + s);
      points->InsertNextPoint(x + s, y, z + s);
      points->InsertNextPoint(x + s, y + s, z + s);
      points->InsertNextPoint(x, y + s, z + s);
    }

    grid->SetPoints(points);

    for (vtkIdType i = 0; i < static_cast<vtkIdType>(mesh.voxels.size()); ++i) {

      vtkIdType ids[8];

      for (int j = 0; j < 8; ++j) {
        ids[j] = 8 * i + j;
      }

      grid->InsertNextCell(VTK_HEXAHEDRON, 8, ids);
    }
  }

public:
  explicit VoxelMeshSaver(const VoxelMesh &mesh,
                          const std::string &path = "voxels.vtu")
      : mesh(mesh), path(path) {

    buildGrid();
  }

  template <typename T>
  void addField(const std::string &name, const std::vector<T> &values) {

    const size_t expectedSize = mesh.voxels.size();

    if (values.size() != expectedSize) {
      throw std::invalid_argument(
          "Field size does not match the number of voxels.");
    }

    addFieldToArray(grid->GetCellData(), name, values);
  }

  void saveMesh() const {

    auto writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();

    writer->SetFileName(path.c_str());
    writer->SetInputData(grid);
    writer->SetDataModeToBinary();

    if (!writer->Write()) {
      throw std::runtime_error("Failed to write voxel mesh to: " + path);
    }
  }
};

#endif
