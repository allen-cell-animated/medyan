
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v3.1
//
//  Copyright (2017-2018)  Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the MEDYAN web page for more information:
//  http://www.medyan.org
//------------------------------------------------------------------

#ifndef MEDYAN_TriangleCylinderBeadExclVolRepulsion_h
#define MEDYAN_TriangleCylinderBeadExclVolRepulsion_h

#include "MathFunctions.h"

//FORWARD DECLARATIONS
class Vertex;
class Bead;

/// Represents a repulsive excluded volume potential used by the
/// TriangleCylinderExclVolume template.
class TriangleCylinderBeadExclVolRepulsion {
    
public:
    double energy(mathfunc::Vec3, mathfunc::Vec3, mathfunc::Vec3, mathfunc::Vec3 cb, double area, double kExVol);
    
    void forces(
        double* f0, double* f1, double* f2, double* fb,
        mathfunc::Vec3, mathfunc::Vec3, mathfunc::Vec3, mathfunc::Vec3 cb,
        double area, const mathfunc::Vec3&, const mathfunc::Vec3&, const mathfunc::Vec3&,
        double kExVol
    );

    mathfunc::Vec3 loadForces(
        const mathfunc::Vec3& c0, const mathfunc::Vec3& c1, const mathfunc::Vec3& c2, const mathfunc::Vec3& coord,
        double area, double kExVol
    );
};


#endif
