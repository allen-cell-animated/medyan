//
//  NeighborList.cpp
//  Cyto
//
//  Created by James Komianos on 11/12/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#include "NeighborList.h"

#include "BeadDB.h"
#include "GController.h"
#include "Cylinder.h"
#include "BoundaryElement.h"

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
                 else if(abs(cylinder->getPositionFilament()
                             - nearbyCylinder->getPositionFilament()) <= 2) continue;
            }
            
            ///Dont add if not within range
            double dist = TwoPointDistance(cylinder->coordinate, nearbyCylinder->coordinate);
            if(dist > _rMax || dist < _rMin) continue;
            
            ///If we got through all of this, add it!
            _list[cylinder].push_back(nearbyCylinder);
        }
    }
}

void CylinderNeighborList::addNeighbor(Neighbor* n) {
    
    ///return if not a cylinder!
    if(!dynamic_cast<Cylinder*>(n)) return;
    ///update neighbors
    updateNeighbors(n);
    
}

void BoundaryElementNeighborList::updateNeighbors(Neighbor* n) {
    
    ///clear existing
    BoundaryElement* be = static_cast<BoundaryElement*>(n);
    _list[be].clear();
    
    ///loop through beads, add as
    for (auto &b : *BeadDB::instance()) {
        
        double dist = be->distance(b->coordinate);
        if(dist > _rMax || dist < _rMin) continue;
        
        ///If we got through this, add it!
        _list[be].push_back(b);
    }
}


void BoundaryElementNeighborList::addNeighbor(Neighbor* n) {
    
    ///return if not a cylinder!
    if(!dynamic_cast<BoundaryElement*>(n)) return;
    ///update neighbors
    updateNeighbors(n);
}





