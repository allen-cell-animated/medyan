
//------------------------------------------------------------------
//  **M3SYM** - Simulation Package for the Mechanochemical
//              Dynamics of Active Networks, 3rd Generation
//
//  Copyright (2014) Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the Papoian lab page for more information:
//  http://papoian.chem.umd.edu/
//------------------------------------------------------------------

#include "NeighborList.h"

#include "Bead.h"
#include "Cylinder.h"
#include "BoundaryElement.h"

#include "GController.h"
#include "MathFunctions.h"

using namespace mathfunc;

void CCNeighborList::addNeighbor(Neighbor* n) {
    
    //return if not a cylinder!
    if(!dynamic_cast<Cylinder*>(n)) return;
    
    //update neighbors
    updateNeighbors(n);
}

void CCNeighborList::updateNeighbors(Neighbor* n) {

    //clear existing
    _list[n].clear();
    Cylinder* cylinder = (Cylinder*)(n);
    
    //Find surrounding compartments (For now its conservative, change soon)
    vector<Compartment*> compartments;
    
    GController::findCompartments(cylinder->coordinate, cylinder->getCompartment(),
            SysParams::Geometry().largestCompartmentSide * 2, compartments);
    
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
            double dist = twoPointDistance(cylinder->coordinate,
                                           nearbyCylinder->coordinate);
            if(dist > _rMax || dist < _rMin) continue;
            
            //If we got through all of this, add it!
            _list[cylinder].push_back(nearbyCylinder);
        }
    }
}

void CCNeighborList::removeDynamicNeighbor(Neighbor* n) {
    
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


vector<Cylinder*> CCNeighborList::getNeighbors(Cylinder* cylinder) {
    
    auto neighbors = _list[cylinder];
    vector<Cylinder*> cylinderNeighbors(_list[cylinder].size());
    
    transform(neighbors.begin(), neighbors.end(), cylinderNeighbors.begin(),
                                   [](Neighbor* n){return (Cylinder*)(n);});
    return cylinderNeighbors;
}

void BBENeighborList::addNeighbor(Neighbor* n) {
    
    //return if not a boundary element!
    if(!dynamic_cast<BoundaryElement*>(n)) return;
    
    //update neighbors
    updateNeighbors(n);
}

void BBENeighborList::updateNeighbors(Neighbor* n) {
    
    //clear existing
    _list[n].clear();
    
    //loop through beads, add as neighbor
    for (auto &b : *BeadDB::instance()) {
        
        double dist = ((BoundaryElement*)(n))->distance(b->coordinate);
        //If within range, add it
        if(dist < _rMax) _list[n].push_back(b);
    }
}

void BBENeighborList::addDynamicNeighbor(Neighbor* n) {
    
    //return if not a boundary element!
    Bead* b; if(!(b = dynamic_cast<Bead*>(n))) return;

    for(auto it = _list.begin(); it != _list.end(); it++) {
        BoundaryElement* be = (BoundaryElement*)(it->first);
        if(be->distance(b->coordinate) < _rMax)
            it->second.push_back(b);
    }
}

void BBENeighborList::removeDynamicNeighbor(Neighbor* n) {
    
    //return if not a boundary element!
    if(!dynamic_cast<Bead*>(n)) return;
    
    for(auto it = _list.begin(); it != _list.end(); it++) {
        auto bead = find(it->second.begin(), it->second.end(), n);
        if(bead != it->second.end()) it->second.erase(bead);
    }
}

vector<Bead*> BBENeighborList::getNeighbors(BoundaryElement* be) {
    
    auto neighbors = _list[be];
    vector<Bead*> beadNeighbors(_list[be].size());
    
    transform(neighbors.begin(), neighbors.end(), beadNeighbors.begin(),
                                   [](Neighbor* n){return (Bead*)(n);});
    return beadNeighbors;
}

void CBENeighborList::addNeighbor(Neighbor* n) {
    
    //return if not a boundary element!
    if(!dynamic_cast<BoundaryElement*>(n)) return;
    
    //update neighbors
    updateNeighbors(n);
}

void CBENeighborList::updateNeighbors(Neighbor* n) {
    
    //clear existing
    _list[n].clear();
    
    //loop through beads, add as neighbor
    for (auto &c : *CylinderDB::instance()) {
        
        double dist = ((BoundaryElement*)(n))->distance(c->coordinate);
        //If within range, add it
        if(dist < _rMax) _list[n].push_back(c);
    }
}

void CBENeighborList::addDynamicNeighbor(Neighbor* n) {
    
    //return if not a boundary element!
    Cylinder* c; if(!(c = dynamic_cast<Cylinder*>(n))) return;
    
    for(auto it = _list.begin(); it != _list.end(); it++) {
        BoundaryElement* be = (BoundaryElement*)(it->first);
        if(be->distance(c->coordinate) < _rMax)
            it->second.push_back(c);
    }
}

void CBENeighborList::removeDynamicNeighbor(Neighbor* n) {
    
    //return if not a boundary element!
    if(!dynamic_cast<Cylinder*>(n)) return;
    
    for(auto it = _list.begin(); it != _list.end(); it++) {
        auto cylinder = find(it->second.begin(), it->second.end(), n);
        if(cylinder != it->second.end()) it->second.erase(cylinder);
    }
}

vector<Cylinder*> CBENeighborList::getNeighbors(BoundaryElement* be) {
    
    auto neighbors = _list[be];
    vector<Cylinder*> cylinderNeighbors(_list[be].size());
    
    transform(neighbors.begin(), neighbors.end(), cylinderNeighbors.begin(),
              [](Neighbor* n){return (Cylinder*)(n);});
    return cylinderNeighbors;
}

