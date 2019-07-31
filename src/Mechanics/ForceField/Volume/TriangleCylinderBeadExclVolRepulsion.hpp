
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

#ifndef MEDYAN_Mechanics_ForceField_Volume_TriangleCylinderBeadExclVolRepulsion_Hpp
#define MEDYAN_Mechanics_ForceField_Volume_TriangleCylinderBeadExclVolRepulsion_Hpp

#include "MathFunctions.h"

/// Represents a repulsive excluded volume potential used by the
/// TriangleBeadExclVolume template.
struct TriangleCylinderBeadExclVolRepulsion {
    
    double energy(mathfunc::Vec3, mathfunc::Vec3, mathfunc::Vec3, mathfunc::Vec3 cb, double area, double kExVol) const;
    
    void forces(
        floatingpoint* f0, floatingpoint* f1, floatingpoint* f2, floatingpoint* fb,
        mathfunc::Vec3, mathfunc::Vec3, mathfunc::Vec3, mathfunc::Vec3 cb,
        double area, const mathfunc::Vec3&, const mathfunc::Vec3&, const mathfunc::Vec3&,
        double kExVol
    ) const;

    mathfunc::Vec3 loadForces(
        const mathfunc::Vec3& c0, const mathfunc::Vec3& c1, const mathfunc::Vec3& c2, const mathfunc::Vec< 3, floatingpoint >& coord,
        double area, double kExVol
    ) const;
};


#endif
