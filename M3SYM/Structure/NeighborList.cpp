
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
    
    ///return if not a cylinder!
    if(!dynamic_cast<Cylinder*>(n)) return;
    
    ///update neighbors
    updateNeighbors(n);
    
}

void CylinderNeighborList::updateNeighbors(Neighbor* n) {

    ///clear existing
    _list[n].clear();
    Cylinder* cylinder = static_cast<Cylinder*>(n);
    
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
            double dist = twoPointDistance(cylinder->coordinate, nearbyCylinder->coordinate);
            if(dist > _rMax || dist < _rMin) continue;
            
            ///If we got through all of this, add it!
            _list[cylinder].push_back(nearbyCylinder);
        }
    }
}

vector<Cylinder*> CylinderNeighborList::getNeighbors(Cylinder* cylinder) {
    
    auto neighbors = _list[cylinder];
    vector<Cylinder*> cylinderNeighbors(_list[cylinder].size());
    
    transform(neighbors.begin(), neighbors.end(), cylinderNeighbors.begin(),
              [](Neighbor* n){return static_cast<Cylinder*>(n);});
    return vector<Cylinder*>(cylinderNeighbors.begin(), cylinderNeighbors.end());
}

void BoundaryElementNeighborList::addNeighbor(Neighbor* n) {
    
    ///return if not a boundary element!
    if(!dynamic_cast<BoundaryElement*>(n)) return;
    
    ///update neighbors
    updateNeighbors(n);
}

void BoundaryElementNeighborList::updateNeighbors(Neighbor* n) {
    
    ///clear existing
    _list[n].clear();
    
    ///loop through beads, add as
    for (auto &b : *BeadDB::instance()) {
        
        double dist = static_cast<BoundaryElement*>(n)->distance(b->coordinate);
        if(dist > _rMax || dist < _rMin) continue;
        
        ///If we got through this, add it!
        _list[n].push_back(b);
    }
}

vector<Bead*> BoundaryElementNeighborList::getNeighbors(BoundaryElement* be) {
    
    auto neighbors = _list[be];
    vector<Bead*> beadNeighbors(_list[be].size());
    
    transform(neighbors.begin(), neighbors.end(), beadNeighbors.begin(),
              [](Neighbor* n){return static_cast<Bead*>(n);});
    return vector<Bead*>(beadNeighbors.begin(), beadNeighbors.end());
}




