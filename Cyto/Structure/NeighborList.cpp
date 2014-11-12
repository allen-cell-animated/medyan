//
//  NeighborList.cpp
//  Cyto
//
//  Created by James Komianos on 11/12/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#include "NeighborList.h"

#include "CylinderDB.h"
#include "GController.h"

#include "MathFunctions.h"

using namespace mathfunc;

vector<Cylinder*> CylinderNeighborList::findNearbyCylinders(Cylinder* cylinder) {
    
    vector<Cylinder*> cylinders;
    
    ///Find surrounding compartments (For now its conservative, change soon)
    vector<Compartment*> compartments;
    GController::findCompartments(cylinder->coordinate, cylinder->getCompartment(),
                 SystemParameters::Geometry().largestCompartmentSide * 2, compartments);
    
    for(auto &c : compartments) {
        for(auto &nearbyCylinder : c->getCylinders()) {
            
            ///dont add cylinders within 1 of this cylinder, on the same filament
            if(cylinder->getFilament() == nearbyCylinder->getFilament()) {
                
                if(_crossFilamentOnly) continue;
                else if(abs(cylinder->getPositionFilament() - nearbyCylinder->getPositionFilament()) <= 1) continue;
            }
            cylinders.push_back(nearbyCylinder);
            
        }
    }
    return vector<Cylinder*>(cylinders.begin(), cylinders.end());
}


void CylinderNeighborList::reset() {
    
    _list.clear(); /// clear current list
    
    //loop through all cylinders
    for(auto &cylinder : *CylinderDB::instance(CylinderDBKey())) {
        
        for(auto &nearbyCylinder : findNearbyCylinders(cylinder)) {
            
            if(cylinder->getID() > nearbyCylinder->getID()) {
                
                double dist = TwoPointDistance(cylinder->coordinate, nearbyCylinder->coordinate);
                if(dist < _rMax && dist > _rMin) _list[cylinder].push_back(nearbyCylinder);
            }
        }
    }
}

void CylinderNeighborList::addNeighbor(Neighbor* n) {
    
    Cylinder* cylinder;
    if((cylinder = dynamic_cast<Cylinder*>(n))) {
        
        for(auto &nearbyCylinder : findNearbyCylinders(cylinder)) {
            
            if(cylinder->getID() > nearbyCylinder->getID()) {
                
                double dist = TwoPointDistance(cylinder->coordinate, nearbyCylinder->coordinate);
                if(dist < _rMax && dist > _rMin) _list[cylinder].push_back(n);
            }
        }
    }
    
}
void CylinderNeighborList::removeNeighbor(Neighbor* n) {
    
    Cylinder* cylinder;
    if((cylinder = dynamic_cast<Cylinder*>(n))) _list.erase(n);
}


