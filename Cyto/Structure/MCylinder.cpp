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
    
    ///set excluded volume const
    SetExVolConst(SystemParameters::Mechanics().VolumeK);
}


void MCylinder::updateExVolNeighborsList(std::vector<MCylinder*>& nearbyMCylinders) {
    
    
    ///Loop through current neighbors, remove any that are not within cutoff
    std::vector<MCylinder*> toRemove;
    for(auto &m : _exVolNeighborsList)
        if(TwoPointDistance(m->_coordinate, _coordinate) >= SystemParameters::Mechanics().VolumeCutoff)
            toRemove.push_back(m);
    
    for(auto &m : toRemove) _exVolNeighborsList.erase(std::find(_exVolNeighborsList.begin(), _exVolNeighborsList.end(), m));
    
    ///Add any that are within cutoff
    for(auto &m : nearbyMCylinders) {
        if(TwoPointDistance(m->_coordinate, _coordinate) < SystemParameters::Mechanics().VolumeCutoff)
            _exVolNeighborsList.push_back(m);
    }
    
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

void MCylinder::SetExVolConst(double k) {_kExVol = k;}
double MCylinder::GetExVolConst() {return _kExVol;}
