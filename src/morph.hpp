#include "mesh_loader.hpp"
#include <vtkUnstructuredGrid.h>

void saveMorphAnimation(const TetrahedralMesh &source,
                        const TetrahedralMesh &target,
                        const std::string &outputName, int nbFrames = 100);
