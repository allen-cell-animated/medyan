#ifndef MEDYAN_MeshSlicing_h
#define MEDYAN_MeshSlicing_h

#include "common.h"

#include <unordered_set>

#include "SimpleGeometricObject.h"

#include "Triangle.h"
#include "GTriangle.h"

/******************************************************************************
Provides classes and functions for slicing through the meshwork.

TODO:
    Finish the implementation, test it and use it
        - OR -
    Delete all of these (including SimpleGeometricObjects.h and stuff)
******************************************************************************/

// Forward declarations
class Membrane;
class Triangle;

class PlaneSliceSnapshot {
private:
    size_t _aspect; // 0: yOz, 1: zOx, 2: xOy
    double _otherCoord; // The remaining coordinate (location) of this snapshot

public:
    PlaneSliceSnapshot(size_t aspect, double otherCoord): _aspect(aspect), _otherCoord(otherCoord) {}

    vector<SimplePolygon<2>> polygons;

    // Properties
    size_t getAspect()const { return _aspect; }

    double getOtherCoord()const { return _otherCoord; }
};


/******************************************************************************
The manager deals with the administrating work associated with mesh slicing.

Provides
    - the simple vertex container.
    - various functions for dealing with slicing
******************************************************************************/
class MeshSlicingManager {
private:
    unordered_set<unique_ptr<SimpleVertex<3>>> vertices3;
    unordered_set<unique_ptr<SimpleVertex<2>>> vertices2;

    // The plane normal from different aspects.
    static const array<array<double, 3>, 3> planeNormal;

    void planeSliceTriangle(size_t aspect, double otherCoord, Triangle* triangle);
    void restoreSlicedTriangle(Triangle* triangle) {
        triangle->getGTriangle()->getPolygons().clear();
        // Vertices constructed would still be hanging around until the manager dies.
    }

    PlaneSliceSnapshot planeSliceMembrane(size_t aspect, double otherCoord, unordered_set<Triangle*>& triangles);
};

#endif
