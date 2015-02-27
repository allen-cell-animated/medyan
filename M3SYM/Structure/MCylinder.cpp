
//------------------------------------------------------------------
//  **M3SYM** - Simulation Package for the Mechanochemical
//              Dynamics of Active Networks, 3rd Generation
//
//  Copyright (2014) Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the Papoian lab page for more information:
//  http://papoian.chem.umd.edu/
//------------------------------------------------------------------

#include "MCylinder.h"

#include "MathFunctions.h"
#include "SysParams.h"

using namespace mathfunc;

MCylinder::MCylinder(double eqLength){
    
    //Set equilibrium length relative to full cylinder length
    setEqLength(eqLength);
    
    //set excluded volume const
    _kExVol = SysParams::Mechanics().VolumeK;
    
    //set angle
    _eqTheta = SysParams::Mechanics().FBendingTheta;
    _eqPhi = SysParams::Mechanics().FTwistingPhi;
}

void MCylinder::setEqLength(double l) {
    _eqLength = l;
    double fracCylinderSize = SysParams::Geometry().cylinderSize / l;
    
    // recalculate other constants
    _kStretch = SysParams::Mechanics().FStretchingK * fracCylinderSize;
    _kBend = SysParams::Mechanics().FBendingK / fracCylinderSize;
}
