#include "common.h"

#include "SimpleGeometricObject.h"

class CompartmentSliceSnapshot {
private:
    size_t _aspect; // 0: yOz, 1: zOx, 2: xOy
    double _otherCoord; // The remaining coordinate (location) of this snapshot

public:

    vector<SimplePolygon<2>> polygons;

    // Properties
    size_t getAspect()const { return _aspect; }

    double getOtherCoord()const { return _otherCoord; }
};