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
    vector<Cylinder*> nearbyCylinders; vector<Compartment*> compartments;
    
    GController::findCompartments(cylinder->coordinate, cylinder->getCompartment(),
            SystemParameters::Geometry().largestCompartmentSide * 2, compartments);
    
    for(auto &c : compartments) {
        for(auto &nearbyCylinder : c->getCylinders()) {
            
            ///dont add cylinders within 1 of this cylinder, on the same filament
            if(cylinder->getFilament() == nearbyCylinder->getFilament()) {
                
                if(_crossFilamentOnly) continue;
                else if(abs(cylinder->getPositionFilament()
                        - nearbyCylinder->getPositionFilament()) <= 1) continue; 
            }
            nearbyCylinders.push_back(nearbyCylinder);
        }
    }

    ///loop through nearby cylinders, add if needed
    for(auto &nearbyCylinder : nearbyCylinders) {
        
        if(cylinder->getID() > nearbyCylinder->getID()) {
            
            double dist = TwoPointDistance(cylinder->coordinate, nearbyCylinder->coordinate);
            if(dist < _rMax && dist > _rMin) _list[cylinder].push_back(nearbyCylinder);
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
    updateNeighbors(n);
    
}
void CylinderNeighborList::removeNeighbor(Neighbor* n) {

    ///return if not a cylinder!
    if(!dynamic_cast<Cylinder*>(n)) return;
    _list.erase(n);
}

const vector<Neighbor*>& CylinderNeighborList::getNeighbors(Neighbor* n) {
    
    return _list[n];
}

