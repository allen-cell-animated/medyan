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
    SetEqLength(eqLength);
    
    _NeighbourList.assign (0, NULL);
}

void MCylinder::SetEqLength(double l) {
    _eqLength = l;
    double fracCylinderSize = SystemParameters::Geometry().cylinderSize / l;
    
    //recalculate other constants
    _kStretch = SystemParameters::Mechanics().FStretchingK * fracCylinderSize;
    _kBend = SystemParameters::Mechanics().FBendingK * fracCylinderSize;
    _kTwist = SystemParameters::Mechanics().FTwistingK * fracCylinderSize;
}

double MCylinder::GetEqLength() {return _eqLength;}

void MCylinder::SetAngle(double alpha) {_eqAngle = alpha;}
double MCylinder::GetAngle() {return _eqAngle;}

void MCylinder::SetStretchingConst(double k) {_kStretch = k;}
double MCylinder::GetStretchingConst() {return _kStretch;}

void MCylinder::SetBendingConst(double k) {_kBend = k;}
double MCylinder::GetBendingConst() {return _kBend;}

void MCylinder::SetTwistingConst(double k) {_kTwist = k;}
double MCylinder::GetTwistingConst() {return _kTwist;}
