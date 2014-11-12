//
//  MCylinder.cpp
//  CytoMech
//
//  Created by Konstantin Popov on 6/30/14.
//  Copyright (c) 2014 Konstantin Popov. All rights reserved.
//
#include <vector>
#include "MathFunctions.h"
#include "MCylinder.h"
#include "Bead.h"
#include "SystemParameters.h"

using namespace mathfunc;

MCylinder::MCylinder(double eqLength){
    
    ///Set equilibrium length relative to full cylinder length
    setEqLength(eqLength);
    
    ///set excluded volume const
    setExVolConst(SystemParameters::Mechanics().VolumeK);
}

void MCylinder::setEqLength(double l) {
    _eqLength = l;
    double fracCylinderSize = SystemParameters::Geometry().cylinderSize / l;
    
    //recalculate other constants
    _kStretch = SystemParameters::Mechanics().FStretchingK * fracCylinderSize;
    _kBend = SystemParameters::Mechanics().FBendingK * fracCylinderSize;
    _kTwist = SystemParameters::Mechanics().FTwistingK * fracCylinderSize;
}
