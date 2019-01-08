
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
    double energy(Vertex*, Vertex*, Vertex*, Bead*, double area, double kExVol, bool stretched);
    
    void forces(
        Vertex*, Vertex*, Vertex*, Bead*,
        double area, const mathfunc::Vec3&, const mathfunc::Vec3&, const mathfunc::Vec3&,
        double kExVol
    );
    void forcesAux(
        Vertex*, Vertex*, Vertex*, Bead*,
        double area, const mathfunc::Vec3&, const mathfunc::Vec3&, const mathfunc::Vec3&,
        double kExVol
    );

    mathfunc::Vec3 loadForces(
        Vertex* v0, Vertex* v1, Vertex* v2, const mathfunc::Vec3& coord,
        double area, double kExVol
    );
};


#endif
