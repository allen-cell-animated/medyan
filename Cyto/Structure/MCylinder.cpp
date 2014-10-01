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

MCylinder::MCylinder(Filament* pf, Bead* firstBead, Bead* secondBead, double eqLength){
    
    ///Set equilibrium length relative to current length
    SetEqLength(eqLength);
    
    _pSecond = secondBead;
    _pFirst = firstBead;
    _NeighbourList.assign (0, NULL);
}

Bead* MCylinder::GetFirstBead() { return _pFirst;}
Bead* MCylinder::GetSecondBead() { return _pSecond;}

void MCylinder::SetEqLength(double l) {
    _eqLength = l;
#ifdef MECHANICS
    //recalculate other constants
    _kStretch = SystemParameters::Mechanics().FStretchingK * SystemParameters::Geometry().cylinderSize / l;
    _kBend = SystemParameters::Mechanics().FBendingK * SystemParameters::Geometry().cylinderSize / l;
    _kTwist = SystemParameters::Mechanics().FTwistingK * SystemParameters::Geometry().cylinderSize / l;
#endif
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
