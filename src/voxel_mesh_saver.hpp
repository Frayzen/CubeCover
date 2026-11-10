#ifndef VOXEL_MESH_SAVER_HPP
#define VOXEL_MESH_SAVER_HPP

#include "ids.hpp"
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

    const size_t nx = mesh.voxels.dimension(0);
    const size_t ny = mesh.voxels.dimension(1);
    const size_t nz = mesh.voxels.dimension(2);

    // Point IDs indexed by grid vertex coordinates.
    Eigen::Tensor<vtkIdType, 3> pointIds(nx + 1, ny + 1, nz + 1);

    pointIds.setConstant(-1);

    // Create each grid point only once.
    auto getPointId = [&](size_t x, size_t y, size_t z) -> vtkIdType {
      vtkIdType &id = pointIds(x, y, z);

      if (id == -1) {
        const double px = mesh.origin.x() + x * s;
        const double py = mesh.origin.y() + y * s;
        const double pz = mesh.origin.z() + z * s;

        id = points->InsertNextPoint(px, py, pz);
      }

      return id;
    };

    // Create hexahedra for occupied voxels.
    for (VoxelId id = 0; id < mesh.nbVoxel; id++) {
      Voxel v = mesh.getVoxel(id);
      auto x = v.x();
      auto y = v.y();
      auto z = v.z();

      if (mesh.voxels(x, y, z) == NO_VOXEL)
        continue;

      vtkIdType ids[8] = {getPointId(x, y, z),
                          getPointId(x + 1, y, z),
                          getPointId(x + 1, y + 1, z),
                          getPointId(x, y + 1, z),

                          getPointId(x, y, z + 1),
                          getPointId(x + 1, y, z + 1),
                          getPointId(x + 1, y + 1, z + 1),
                          getPointId(x, y + 1, z + 1)};

      grid->InsertNextCell(VTK_HEXAHEDRON, 8, ids);
    }

    grid->SetPoints(points);
  }

public:
  explicit VoxelMeshSaver(const VoxelMesh &mesh,
                          const std::string &path = "voxels.vtu")
      : mesh(mesh), path(path) {

    buildGrid();
  }

  template <typename T>
  void addField(const std::string &name, const std::vector<T> &values) {

    const size_t expectedSize = mesh.nbVoxel;

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
