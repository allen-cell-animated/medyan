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
    auto& po = triangle->getGTriangle()->getPolygons();
    if(po.size() == 0) {
        // Create the first polygon using this triangle
        auto sv0 = make_unique<SimpleVertex<3>>(vector2Array<double, 3>(v0));
        auto sv1 = make_unique<SimpleVertex<3>>(vector2Array<double, 3>(v1));
        auto sv2 = make_unique<SimpleVertex<3>>(vector2Array<double, 3>(v2));
        po.emplace_back({{sv0.get(), sv1.get(), sv2.get()}});

        vertices3.push_back(move(sv0));
        vertices3.push_back(move(sv1));
        vertices3.push_back(move(sv2));
    }

    for(size_t idx = 0; idx < po.size(); ++idx) { // Note that now size() could change.

    }
    // TODO: implementation
}

