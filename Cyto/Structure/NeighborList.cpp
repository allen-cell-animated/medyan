//
//  NeighborList.cpp
//  Cyto
//
//  Created by James Komianos on 11/12/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#include "NeighborList.h"

#include "Cylinder.h"
#include "GController.h"

#include "MathFunctions.h"

using namespace mathfunc;

void CylinderNeighborList::updateNeighbors(Neighbor* n) {

    ///clear existing
    Cylinder* cylinder = static_cast<Cylinder*>(n);
    _list[cylinder].clear();
    
    ///Find surrounding compartments (For now its conservative, change soon)
    vector<Compartment*> compartments;
    
    GController::findCompartments(cylinder->coordinate, cylinder->getCompartment(),
            SystemParameters::Geometry().largestCompartmentSide * 2, compartments);
    
    for(auto &c : compartments) {
        for(auto &nearbyCylinder : c->getCylinders()) {
            
            ///Dont add if ID is more than cylinder
            if(cylinder->getID() <= nearbyCylinder->getID()) continue;
            
            ///Don't add if on the same filament
            if(cylinder->getFilament() == nearbyCylinder->getFilament()) {
                
                 ///if cross filament only interaction, dont add
                 if(_crossFilamentOnly) continue;
                
                 ///if not cross filament, check if not neighboring
                 else if(abs(cylinder->getPositionFilament() - nearbyCylinder->getPositionFilament()) <= 2) continue;
            }
            
            ///Dont add if not within range
            double dist = TwoPointDistance(cylinder->coordinate, nearbyCylinder->coordinate);
            if(dist > _rMax || dist < _rMin) continue;
            
            ///If we got through all of this, add it!
            _list[cylinder].push_back(nearbyCylinder);
        }
    }
}

void CylinderNeighborList::reset() {
    
    //loop through all neighbor keys
    for(auto it = _list.begin(); it != _list.end(); it++) {
        
        it->second.clear(); ///clear vector of neighbors
        updateNeighbors(it->first);
    }
}

void CylinderNeighborList::addNeighbor(Neighbor* n) {
    
    ///return if not a cylinder!
    if(!dynamic_cast<Cylinder*>(n)) return;
    ///update neighbors
    updateNeighbors(n);
    
}
void CylinderNeighborList::removeNeighbor(Neighbor* n) {

    ///erase
    _list.erase(n);
}

const vector<Neighbor*>& CylinderNeighborList::getNeighbors(Neighbor* n) {
    
    return _list[n];
}

