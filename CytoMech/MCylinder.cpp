//
//  MCylinder.cpp
//  CytoMech
//
//  Created by Konstantin Popov on 6/30/14.
//  Copyright (c) 2014 Konstantin Popov. All rights reserved.
//
#include <vector>
#include "MCylinder.h"
#include "MBead.h"


MCylinder::MCylinder(Filament* pf, Bead* pb){
    
    _eqLength = SystemParameters::Mechanics().FStretchingL;
    _eqAngle = SystemParameters::Mechanics().FBendingTheta;
    _eqAngleTwist = SystemParameters::Mechanics().FTwistingPhi;
    _kStretch = SystemParameters::Mechanics().FStretchingK;
    _kBend = SystemParameters::Mechanics().FBendingK;
    _kTwist = SystemParameters::Mechanics().FTwistingK;
    
    _pSecond = NULL;
    
    _pFirst = pb;
    _NeighbourList.assign (0, NULL);
}


void MCylinder::SetSecondBead(Bead *pb) {_pSecond = pb;}

Bead* MCylinder::GetFirstBead() { return _pFirst;}
Bead* MCylinder::GetSecondBead() { return _pSecond;}

void MCylinder::SetEqLength(double l) {_eqLength = l;}
double MCylinder::GetEqLength() {return _eqLength;}

void MCylinder::SetAngle(double alpha) {_eqAngle = alpha;}
double MCylinder::GetAngle() {return _eqAngle;}

void MCylinder::SetStretchingConst(double k) {_kStretch = k;}
double MCylinder::GetStretchingConst() {return _kStretch;}

void MCylinder::SetBendingConst(double k) {_kBend = k;}
double MCylinder::GetBendingConst() {return _kBend;}

void MCylinder::SetTwistingConst(double k) {_kTwist = k;}
double MCylinder::GetTwistingConst() {return _kTwist;}
