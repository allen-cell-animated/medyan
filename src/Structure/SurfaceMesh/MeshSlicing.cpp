#include "MeshSlicing.h"

#include "SimpleGeometricObject.h"
#include "Membrane.h"
#include "Triangle.h"
#include "Vertex.h"

#include "MathFunctions.h"
using namespace mathfunc;

void MeshSlicingManager::planeSliceTriangle(size_t aspect, double otherCoord, Triangle* triangle) {
    /**************************************************************************
    After slicing, several things will be produced:
        - Triangles will be sliced and contain 0 or more more 3d polygons
        - A plane slicing snapshot will be created
    **************************************************************************/
    
    // Check whether this triangle can be cut
    auto& v0 = triangle->getVertices()[0]->coordinate;
    auto& v1 = triangle->getVertices()[1]->coordinate;
    auto& v2 = triangle->getVertices()[2]->coordinate;
    if((v0[aspect] - otherCoord) * (v1[aspect] - otherCoord) > 0
        && (v1[aspect] - otherCoord) * (v2[aspect] - otherCoord) > 0)
        return; // No need to cut this triangle
    
    // Check whether this triangle has already been cut
    // TODO: implementation
}

