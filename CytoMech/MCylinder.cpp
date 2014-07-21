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


Cylinder::Cylinder(Filament* pf, Bead* pb){
    
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

bool Cylinder::IfLast(){
    if (_ifLast) return true;
    
    return false;
}

void Cylinder::SetLast(bool b){ _ifLast = b;}

void Cylinder::SetSecondBead(Bead *pb) {_pSecond = pb;}

Bead* Cylinder::GetFirstBead() { return _pFirst;}
Bead* Cylinder::GetSecondBead() { return _pSecond;}

void Cylinder::SetEqLength(double l) {_eqLenght = l;}
double Cylinder::GetEqLength() {return _eqLenght;}

void Cylinder::SetAngle(double alpha) {_eqAngle = alpha;}
double Cylinder::GetAngle() {return _eqAngle;}

void Cylinder::SetStretchingConst(double k) {_kStretch = k;}
double Cylinder::GetStretchingConst() {return _kStretch;}

void Cylinder::SetBendingConst(double k) {_kBend = k;}
double Cylinder::GetBendingConst() {return _kBend;}

void Cylinder::SetTwistingConst(double k) {_kTwist = k;}
double Cylinder::GetTwistingConst() {return _kTwist;}
