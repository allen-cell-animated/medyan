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

void MCylinder::addExVolNeighbor(MCylinder* neighbor) { _exVolNeighborsList.push_back(neighbor);}

void MCylinder::removeExVolNeighbor(MCylinder* neighbor) {
    auto it = find(_exVolNeighborsList.begin(), _exVolNeighborsList.end(), neighbor);
    if(it != _exVolNeighborsList.end()) _exVolNeighborsList.erase(it);
}


MCylinder::MCylinder(double eqLength){
    
    ///Set equilibrium length relative to full cylinder length
    setEqLength(eqLength);
    
    ///set excluded volume const
    setExVolConst(SystemParameters::Mechanics().VolumeK);
}

MCylinder::~MCylinder() {
    ///Remove from current neighbors excluded volume list
    for(auto &neighbor : _exVolNeighborsList) neighbor->removeExVolNeighbor(this);
}

void MCylinder::updateExVolNeighborsList(vector<MCylinder*>& nearbyMCylinders) {
    
    ///Loop through current neighbors, remove any that are not within cutoff
    vector<MCylinder*> toRemove;
    for(auto &m : _exVolNeighborsList)
        if(TwoPointDistance(m->_coordinate, _coordinate) >= SystemParameters::Mechanics().VolumeCutoff) toRemove.push_back(m);
    
    for(auto &m : toRemove) removeExVolNeighbor(m);
    
    ///Add any that are within cutoff
    for(auto &m : nearbyMCylinders)
        if(TwoPointDistance(m->_coordinate, _coordinate) < SystemParameters::Mechanics().VolumeCutoff) addExVolNeighbor(m);
    
}

void MCylinder::setEqLength(double l) {
    _eqLength = l;
    double fracCylinderSize = SystemParameters::Geometry().cylinderSize / l;
    
    //recalculate other constants
    _kStretch = SystemParameters::Mechanics().FStretchingK * fracCylinderSize;
    _kBend = SystemParameters::Mechanics().FBendingK * fracCylinderSize;
    _kTwist = SystemParameters::Mechanics().FTwistingK * fracCylinderSize;
}
