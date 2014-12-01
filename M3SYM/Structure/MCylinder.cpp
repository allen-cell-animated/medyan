
//------------------------------------------------------------------
//  **M3SYM** - Simulation Package for the Mechanochemical
//              Dynamics of Active Networks, 3rd Generation
//
//  Copyright (2014) Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the Papoian lab page for installation and documentation:
//  http://papoian.chem.umd.edu/
//------------------------------------------------------------------

#include "MCylinder.h"

#include "MathFunctions.h"
#include "SystemParameters.h"

using namespace mathfunc;

MCylinder::MCylinder(double eqLength){
    
    //Set equilibrium length relative to full cylinder length
    setEqLength(eqLength);
    
    //set excluded volume const
    setExVolConst(SystemParameters::Mechanics().VolumeK);
}

void MCylinder::setEqLength(double l) {
    _eqLength = l;
    double fracCylinderSize = SystemParameters::Geometry().cylinderSize / l;
    
    // recalculate other constants
    _kStretch = SystemParameters::Mechanics().FStretchingK * fracCylinderSize;
    _kBend = SystemParameters::Mechanics().FBendingK * fracCylinderSize;
    _kTwist = SystemParameters::Mechanics().FTwistingK * fracCylinderSize;
}
