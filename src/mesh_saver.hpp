#ifndef MESH_SAVER_HPP
#define MESH_SAVER_HPP

#include "mesh_loader.hpp"
#include "vertex.hpp"

#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

#include <Eigen/Dense>

#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkUnsignedCharArray.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkXMLUnstructuredGridWriter.h>

enum class ElementType { Node, Edge, Cell };

std::vector<Edge> buildEdges(const std::vector<Tetrahedron> &tets);

class MeshSaver {
private:
  vtkSmartPointer<vtkUnstructuredGrid> grid;
  vtkSmartPointer<vtkPolyData> edgeGrid;

  std::string path;

  size_t nbNodes;
  size_t nbEdges;
  size_t nbCells;

  template <typename T> static constexpr int nbComponents() {
    if constexpr (std::is_same_v<T, bool>) {
      return 1;
    } else if constexpr (std::is_arithmetic_v<T>) {
      return 1;
    } else if constexpr (requires {
                           T::RowsAtCompileTime;
                           T::ColsAtCompileTime;
                         }) {
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
    array->SetNumberOfTuples(values.size());

    for (size_t i = 0; i < values.size(); ++i) {
      double tuple[nComp];

      for (int j = 0; j < nComp; ++j) {
        tuple[j] = component(values[i], j);
      }

      array->SetTuple(i, tuple);
    }

    attributes->AddArray(array);
  }

public:
  MeshSaver(const std::vector<Vertex> &vertices, const std::vector<Edge> &edges,
            const std::vector<Tetrahedron> &tets, const std::string &path)
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

  MeshSaver(const std::vector<Vertex> &vertices,
            const std::vector<Tetrahedron> &tets, const std::string &path)
      : MeshSaver(vertices, buildEdges(tets), tets, path) {}

  template <ElementType E, typename T>
  void addField(const std::string &name, const std::vector<T> &values) {
    vtkDataSetAttributes *attributes = nullptr;
    size_t expectedSize = 0;

    if constexpr (E == ElementType::Node) {
      attributes = grid->GetPointData();
      expectedSize = nbNodes;
    } else if constexpr (E == ElementType::Cell) {
      attributes = grid->GetCellData();
      expectedSize = nbCells;
    } else if constexpr (E == ElementType::Edge) {
      attributes = edgeGrid->GetCellData();
      expectedSize = nbEdges;
    }

    if (values.size() != expectedSize) {
      throw std::invalid_argument(
          "Field size does not match the number of elements.");
    }

    addFieldToArray(attributes, name, values);
  }

  void save() const {
    auto writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();

    writer->SetFileName(path.c_str());
    writer->SetInputData(grid);
    writer->SetDataModeToBinary();

    if (writer->Write() != 1) {
      throw std::runtime_error("Failed to write VTK mesh to " + path);
    }

    if (nbEdges > 0) {
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
  }
};
#endif
