
//------------------------------------------------------------------
//  **M3SYM** - Simulation Package for the Mechanochemical
//              Dynamics of Active Networks, 3rd Generation
//
//  Copyright (2014) Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the Papoian lab page for installation and documentation:
//  http://papoian.chem.umd.edu/
//------------------------------------------------------------------

#include "NeighborList.h"

#include "Bead.h"
#include "Cylinder.h"
#include "BoundaryElement.h"

#include "GController.h"
#include "MathFunctions.h"

using namespace mathfunc;

void CylinderNeighborList::addNeighbor(Neighbor* n) {
    
    //return if not a cylinder!
    if(!dynamic_cast<Cylinder*>(n)) return;
    
    //update neighbors
    updateNeighbors(n);
}

void CylinderNeighborList::updateNeighbors(Neighbor* n) {

    //clear existing
    _list[n].clear();
    Cylinder* cylinder = (Cylinder*)(n);
    
    //Find surrounding compartments (For now its conservative, change soon)
    vector<Compartment*> compartments;
    
    GController::findCompartments(cylinder->coordinate, cylinder->getCompartment(),
            SystemParameters::Geometry().largestCompartmentSide * 2, compartments);
    
    for(auto &c : compartments) {
        for(auto &nearbyCylinder : c->getCylinders()) {
            
            //Dont add if ID is more than cylinder
            if(cylinder->getID() <= nearbyCylinder->getID()) continue;
            
            //Don't add if on the same filament
            if(cylinder->getFilament() == nearbyCylinder->getFilament()) {
                
                 //if cross filament only interaction, dont add
                 if(_crossFilamentOnly) continue;
                
                 //if not cross filament, check if not neighboring
                 else if(abs(cylinder->getPositionFilament() -
                   nearbyCylinder->getPositionFilament()) <= 2) continue;
            }
            
            //Dont add if not within range
            double dist = twoPointDistance(cylinder->coordinate, nearbyCylinder->coordinate);
            if(dist > _rMax || dist < _rMin) continue;
            
            //If we got through all of this, add it!
            _list[cylinder].push_back(nearbyCylinder);
        }
    }
}

void CylinderNeighborList::removeDynamicNeighbor(Neighbor* n) {
    
    //return if not a cylinder!
    if(!dynamic_cast<Cylinder*>(n)) return;
    
    //remove its own list
    removeNeighbor(n);
    
    //remove from other lists
    for(auto it = _list.begin(); it != _list.end(); it++) {
        auto cylinder = find(it->second.begin(), it->second.end(), n);
        if(cylinder != it->second.end()) it->second.erase(cylinder);
    }
}


vector<Cylinder*> CylinderNeighborList::getNeighbors(Cylinder* cylinder) {
    
    auto neighbors = _list[cylinder];
    vector<Cylinder*> cylinderNeighbors(_list[cylinder].size());
    
    transform(neighbors.begin(), neighbors.end(), cylinderNeighbors.begin(),
                                   [](Neighbor* n){return (Cylinder*)(n);});
    return cylinderNeighbors;
}

void BoundaryElementNeighborList::addNeighbor(Neighbor* n) {
    
    //return if not a boundary element!
    if(!dynamic_cast<BoundaryElement*>(n)) return;
    
    //update neighbors
    updateNeighbors(n);
}

void BoundaryElementNeighborList::updateNeighbors(Neighbor* n) {
    
    //clear existing
    _list[n].clear();
    
    //loop through beads, add as neighbor
    for (auto &b : *BeadDB::instance()) {
        
        double dist = ((BoundaryElement*)(n))->distance(b->coordinate);
        //If within range, add it
        if(dist < _rMax) _list[n].push_back(b);
    }
}

void BoundaryElementNeighborList::addDynamicNeighbor(Neighbor* n) {
    
    //return if not a boundary element!
    if(!dynamic_cast<Bead*>(n)) return;
    
    //cast to bead, add in each boundary element
    Bead* b = (Bead*)n;
    
    for(auto it = _list.begin(); it != _list.end(); it++) {
        BoundaryElement* be = (BoundaryElement*)(it->first);
        if(be->distance(b->coordinate) < _rMax)
            it->second.push_back(b);
    }
}

void BoundaryElementNeighborList::removeDynamicNeighbor(Neighbor* n) {
    
    //return if not a boundary element!
    if(!dynamic_cast<Bead*>(n)) return;
    
    for(auto it = _list.begin(); it != _list.end(); it++) {
        auto bead = find(it->second.begin(), it->second.end(), n);
        if(bead != it->second.end()) it->second.erase(bead);
    }
}

vector<Bead*> BoundaryElementNeighborList::getNeighbors(BoundaryElement* be) {
    
    auto neighbors = _list[be];
    vector<Bead*> beadNeighbors(_list[be].size());
    
    transform(neighbors.begin(), neighbors.end(), beadNeighbors.begin(),
                                   [](Neighbor* n){return (Bead*)(n);});
    return beadNeighbors;
}

