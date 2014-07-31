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
    
    _eqLenght = L;
    _eqAngle = teta;
    _eqAngleTwist = 0.0;
    _kStretch = kS;
    _kBend = kB;
    _kTwist = kTw;
    
    _pSecond = NULL;
    _ifLast = true;
    
    _positionFilament = -1000;
    
    _pFilament = pf;
    _pFirst = pb;
    _pFirst->SetParent(this);
    _NeighbourList.assign (0, NULL);
}

bool MCylinder::IfLast(){
    if (_ifLast) return true;
    
    return false;
}

void MCylinder::SetLast(bool b){ _ifLast = b;}

void MCylinder::SetSecondBead(Bead *pb) {_pSecond = pb;}

Bead* MCylinder::GetFirstBead() { return _pFirst;}
Bead* MCylinder::GetSecondBead() { return _pSecond;}

void MCylinder::SetEqLength(double l) {_eqLenght = l;}
double MCylinder::GetEqLength() {return _eqLenght;}

void MCylinder::SetAngle(double alpha) {_eqAngle = alpha;}
double MCylinder::GetAngle() {return _eqAngle;}

void MCylinder::SetStretchingConst(double k) {_kStretch = k;}
double MCylinder::GetStretchingConst() {return _kStretch;}

void MCylinder::SetBendingConst(double k) {_kBend = k;}
double MCylinder::GetBendingConst() {return _kBend;}

void MCylinder::SetTwistingConst(double k) {_kTwist = k;}
double MCylinder::GetTwistingConst() {return _kTwist;}
