#ifndef MEDYAN_MeshSlicing_h
#define MEDYAN_MeshSlicing_h

#include "common.h"

#include "SimpleGeometricObject.h"

/******************************************************************************
Provides classes and functions for slicing through the meshwork.
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
public:
    unordered_set<unique_ptr<SimpleVertex<3>>> vertices3;
    unordered_set<unique_ptr<SimpleVertex<2>>> vertices2;

    void planeSliceTriangle(size_t aspect, double otherCoord, Triangle* triangle);
}

#endif
