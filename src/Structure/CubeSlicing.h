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
    // x = x_min, x = x_max,
    // y = y_min, y = y_max,
    // z = z_min, z = z_max

    PlaneCubeSlicingResult& operator*=(double a) {
        // Expand the volume and area
        double a2 = a * a;
        volumeIn *= a2 * a;
        for(double& eachAreaIn: areaIn) eachAreaIn *= a2;
        return *this;
    }

    PlaneCubeSlicingResult& flip(unsigned short int flippedCoords) { // coords (0-7) in binary is zyx
        for(size_t inspect = 0; inspect < 3; ++inspect) {
            if((flippedCoords >> inspect) & 1) {
                std::swap(areaIn[inspect*2], areaIn[inspect*2+1]);
            }
        }
        return *this;
    }
    PlaneCubeSlicingResult& reverse(double a) { // a is cube size
        double a2 = a * a;
        volumeIn = a2 * a - volumeIn;
        for(double& eachAreaIn: areaIn) eachAreaIn = a2 - eachAreaIn;
        return *this;
    }

};

inline PlaneCubeSlicingResult planeUnitCubeSlice( // Cube [0, 1] x [0, 1] x [0, 1]
    const std::array<double, 3>& point,  // A point on the plane
    const std::array<double, 3>& normal, // Unit normal of the plane pointing outwards
) {
    PlaneCubeSlicingResult res{}; // zero-initiate

    // First consider the case where normal has only non-neg components
    std::array<double, 3> flippedPoint = point;
    std::array<double, 3> flippedNormal = normal;
    unsigned short int flippedCoords = 0;
    for(size_t idx = 0; idx < 3; ++idx) {
        if(flippedNormal[idx] < 0) {
            flippedCoords += (1 << idx);
            flippedPoint[idx] = 1 - flippedPoint[idx];
            flippedNormal[idx] = -flippedNormal[idx];
        }
    }

    double signedDist = dotProduct(flippedPoint, flippedNormal);
    double signedDistMax = dotProduct(std::array<double, 3>{{1, 1, 1}}, flippedNormal);

    bool reverse = false;

    if(signedDist <= 0) {
        // do nothing because everything is zero
    } else if(signedDist >= signedDistMax) {
        res.volumeIn = 1;
        res.areaIn = {{1, 1, 1, 1, 1, 1}};
    } else {
        if(signedDist > 0.5 * signedDistMax) {
            // Doing a full flip, plus swapping inside/outside,
            // which leaves the normal unchanged.
            reverse = true;
            flippedPoint ^= ((1 << 3) - 1); // full flip
            for(size_t idx = 0; idx < 3; ++idx) {
                flippedPoint[idx] = 1 - flippedPoint[idx];
            }

            signedDist = signedDistMax - signedDist;
        }

        // Now the signedDist should be within [0, 0.5 signedDistMax]

        // Find the intersections on the axis. The values should be within [0, +inf]
        double x = signedDistance / flippedNormal[0];
        double y = signedDistance / flippedNormal[1];
        double z = signedDistance / flippedNormal[2];

        size_t numGT1 = 0;
        if(x > 1) ++numGT1;
        if(y > 1) ++numGT1;
        if(z > 1) ++numGT1;

        switch(numGT1) {
        case 0:
            res.volumeIn = x * y * z / 6;
            res.areaIn[0] = y * z / 2;
            res.areaIn[2] = z * x / 2;
            res.areaIn[4] = x * y / 2;
            break;

        case 1:
            if(x > 1) {
            } else if(y > 1) {
            } else if(z > 1) {
            }
            break;

        case 2:
            if(x <= 1) {
            } else if(y <= 1) {
            } else if(z <= 1) {
            }
            break;
        
        case 3:
            break;
            // TODO: implementation
        }
    }

    // Finally restore the flipping
    res.flip(flippedCoords);
    if(reverse) res.reverse(1);

    return res;
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
