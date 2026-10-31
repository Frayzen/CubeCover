#ifndef MESH_SAVER_HPP
#define MESH_SAVER_HPP

#include "mesh_loader.hpp"

#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

#include <Eigen/Dense>

#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkDataSetAttributes.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkXMLUnstructuredGridWriter.h>

enum class ElementType { NodeType, EdgeType, FaceType, CellType };

class MeshSaver {
private:
  const TetrahedralMesh &mesh;

  vtkSmartPointer<vtkUnstructuredGrid> grid;
  vtkSmartPointer<vtkPolyData> edgeGrid;
  vtkSmartPointer<vtkPolyData> faceGrid;

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

      const T value = values[i];

      double tuple[nComp];

      for (int j = 0; j < nComp; ++j) {
        tuple[j] = component(value, j);
      }

      array->SetTuple(static_cast<vtkIdType>(i), tuple);
    }

    attributes->AddArray(array);
  }

  void buildVolumeGrid();
  void buildEdgeGrid();
  void buildFaceGrid();

public:
  explicit MeshSaver(const TetrahedralMesh &mesh, const std::string &path);

  template <ElementType E, typename T>
  void addField(const std::string &name, const std::vector<T> &values) {

    vtkDataSetAttributes *attributes = nullptr;
    size_t expectedSize = 0;

    if constexpr (E == ElementType::NodeType) {

      attributes = grid->GetPointData();
      expectedSize = mesh.vertices.size();

    } else if constexpr (E == ElementType::EdgeType) {

      attributes = edgeGrid->GetCellData();
      expectedSize = mesh.nbBoundaryEdge;

    } else if constexpr (E == ElementType::FaceType) {

      attributes = faceGrid->GetCellData();
      expectedSize = mesh.nbBoundaryFace;

    } else if constexpr (E == ElementType::CellType) {

      attributes = grid->GetCellData();
      expectedSize = mesh.tets.size();

    } else {

      static_assert(E == ElementType::NodeType || E == ElementType::EdgeType ||
                        E == ElementType::FaceType ||
                        E == ElementType::CellType,
                    "Unsupported element type");
    }

    if (values.size() != expectedSize) {
      throw std::invalid_argument(
          "Field size does not match the number of elements.");
    }

    addFieldToArray(attributes, name, values);
  }

  void saveMesh() const;
  void saveEdges() const;
  void saveFaces() const;
};
#endif
