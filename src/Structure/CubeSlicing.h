#ifndef MEDYAN_CubeSlicing_h
#define MEDYAN_CubeSlicing_h

#include <array>

#include "MathFunctions.h"
using namespace mathfunc;

//FORWARD DECLARATIONS

// The result produced by plane intersecting a cube
struct PlaneCubeSlicingResult {
    double volumeIn;

    std::array<double, 6> areaIn;
    // in the order of
    // (x_min, y, z), (x_max, y, z),
    // (x, y_min, z), (x, y_max, z),
    // (x, y, z_min), (x, y, z_max)

    PlaneCubeSlicingResult& operator*=(double a) {
        // Expand the volume and area
        double a2 = a * a;
        volumeIn *= a2 * a;
        for(double& eachAreaIn: areaIn) eachAreaIn *= a2;
        return *this;
    }

};

inline PlaneCubeSlicingResult planeUnitCubeSlice( // Cube [0, 1] x [0, 1] x [0, 1]
    const std::array<double, 3>& point,  // A point on the plane
    const std::array<double, 3>& normal, // Unit normal of the plane pointing outwards
) {
    PlaneCubeSlicingResult res;
    // TODO: implementation
}

inline PlaneCubeSlicingResult planeCubeSlice(
    const std::array<double, 3>& point,  // A point on the plane
    const std::array<double, 3>& normal, // Unit normal of the plane pointing outwards
    const std::array<double, 3>& r0,     // (x_min, y_min, z_min) of the cube
    double                       a       // Side length of the cube
) {
    PlaneCubeSlicingResult res;
    res = planeUnitCubeSlice(
        vectorMultiply(vectorDifference(point, r0), 1.0 / a),
        normal
    );
    res *= a;
    return res;
}

#endif
