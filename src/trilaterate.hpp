#ifndef TRILATERATE_HPP
#define TRILATERATE_HPP

#include "vertex.hpp"

// We project the oldTip to the newFace
Vertex tipProject(const Vertex oldFace[3], const Vertex newFace[3],
                  Vertex oldTip);
#endif
