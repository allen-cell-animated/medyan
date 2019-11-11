
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v3.1
//
//  Copyright (2015-2016)  Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the MEDYAN web page for more information:
//  http://www.medyan.org
//------------------------------------------------------------------

#include "NeighborListImpl.h"

#include "Bead.h"
#include "Filament.h"
#include "Cylinder.h"
#include "Bubble.h"
#include "BoundaryElement.h"
#ifdef CAMKII_ENABLED
#include "CaMKIIingPoint.h"
#endif

#include "GController.h"
#include "MathFunctions.h"

using namespace mathfunc;


//CYLINDER-CYLINDER

void CylinderCylinderNL::updateNeighbors(Cylinder* cylinder, bool runtime) {
    
    //clear existing
    _list[cylinder].clear();
    
    //Find surrounding compartments (For now its conservative)
    vector<Compartment*> compartments;
    auto searchDist = SysParams::Geometry().largestCompartmentSide;
    
    GController::findCompartments(cylinder->coordinate,
                                  cylinder->getCompartment(),
                                  searchDist + _rMax, compartments);
    
    for(auto &comp : compartments) {
        for(auto &ncylinder : comp->getCylinders()) {

      //Don't add the same cylinder!
      if(cylinder == ncylinder) continue;

			//Dont add if ID is more than cylinder for half-list
			if (!_full && cylinder->getID() <= ncylinder->getID()) continue;

            //Don't add if belonging to same parent
            if(cylinder->getParent() == ncylinder->getParent()) {

                //if not cross filament, check if not neighboring
                auto dist = fabs(cylinder->getPosition() -
                                 ncylinder->getPosition());
                if(dist <= 2) continue;
            }

            //Dont add if not within range
            double dist = twoPointDistance(cylinder->coordinate,
                                           ncylinder->coordinate);
            if(dist > _rMax || dist < _rMin) continue;

				//If we got through all of this, add it!
				_list[cylinder].push_back(ncylinder);

				//if runtime, add to other list as well if full
				if ((runtime && _full))
					_list[ncylinder].push_back(cylinder);
//			}
        }
    }
}

void CylinderCylinderNL::addNeighbor(Neighbor* n) {

    //return if not a cylinder!
    Cylinder* cylinder;
    if(!(cylinder = dynamic_cast<Cylinder*>(n))) return;
    
    //update neighbors
    updateNeighbors(cylinder, true);
}

void CylinderCylinderNL::removeNeighbor(Neighbor* n) {
    
    Cylinder* cylinder;
    if(!(cylinder = dynamic_cast<Cylinder*>(n))) return;
    
    _list.erase(cylinder);
    
    //remove from other lists
    for(auto it = _list.begin(); it != _list.end(); it++) {
        
        auto cit = find(it->second.begin(), it->second.end(), cylinder);
        if(cit != it->second.end()) it->second.erase(cit);
    }
}

void CylinderCylinderNL::reset() {
    
    _list.clear();
    
    //loop through all neighbor keys
    for(auto cylinder: Cylinder::getCylinders())
        
        updateNeighbors(cylinder);
}

vector<Cylinder*> CylinderCylinderNL::getNeighbors(Cylinder* cylinder) {
    return _list[cylinder];
}

//BOUNDARYELEMENT - CYLINDER

void BoundaryCylinderNL::updateNeighbors(BoundaryElement* be) {
    
    //clear existing
    _list[be].clear();
    
    //loop through beads, add as neighbor
    for (auto &c : Cylinder::getCylinders()) {
        
        double dist = be->distance(c->coordinate);
        //If within range, add it
        if(dist < _rMax) _list[be].push_back(c);
    }
}

void BoundaryCylinderNL::addNeighbor(Neighbor* n) {
    
    //return if not a boundary element!
    BoundaryElement* be;
    if(!(be = dynamic_cast<BoundaryElement*>(n))) return;
    
    //update neighbors
    updateNeighbors(be);
}

void BoundaryCylinderNL::removeNeighbor(Neighbor* n) {
    
    BoundaryElement* be;
    if(!(be = dynamic_cast<BoundaryElement*>(n))) return;
    
    _list.erase(be);
}

void BoundaryCylinderNL::addDynamicNeighbor(DynamicNeighbor* n) {
    
    //return if not a cylinder!
    Cylinder* c;
    
    if(!(c = dynamic_cast<Cylinder*>(n))) return;
    
    for(auto it = _list.begin(); it != _list.end(); it++) {
        
        //if within range, add it
        if(it->first->distance(c->coordinate) < _rMax)
            it->second.push_back(c);
    }
}

void BoundaryCylinderNL::removeDynamicNeighbor(DynamicNeighbor* n) {
    
    //return if not a cylinder!
    Cylinder* c;
    
    if(!(c = dynamic_cast<Cylinder*>(n))) return;
    
    for(auto it = _list.begin(); it != _list.end(); it++) {
        
        auto cit = find(it->second.begin(), it->second.end(), c);
        if(cit != it->second.end()) it->second.erase(cit);
    }
}

void BoundaryCylinderNL::reset() {
    
    _list.clear();
    
    //loop through all neighbor keys
    for(auto boundary: BoundaryElement::getBoundaryElements())
        
        updateNeighbors(boundary);
}

vector<Cylinder*> BoundaryCylinderNL::getNeighbors(BoundaryElement* be) {
    
    return _list[be];
}

//BOUNDARYELEMENT - BUBBLE

void BoundaryBubbleNL::updateNeighbors(BoundaryElement* be) {
    
    //clear existing
    _list[be].clear();
    
    //loop through beads, add as neighbor
    for (auto &b : Bubble::getBubbles()) {
        
        double dist = be->distance(b->coordinate);
        //If within range, add it
        if(dist < _rMax) _list[be].push_back(b);
    }
}

void BoundaryBubbleNL::addNeighbor(Neighbor* n) {
    
    //return if not a boundary element!
    BoundaryElement* be;
    if(!(be = dynamic_cast<BoundaryElement*>(n))) return;
    
    //update neighbors
    updateNeighbors(be);
}

void BoundaryBubbleNL::removeNeighbor(Neighbor* n) {
    
    BoundaryElement* be;
    if(!(be = dynamic_cast<BoundaryElement*>(n))) return;
    
    _list.erase(be);
}

void BoundaryBubbleNL::addDynamicNeighbor(DynamicNeighbor* n) {
    
    //return if not a filament bead!
    Bubble* b;
    
    if(!(b = dynamic_cast<Bubble*>(n))) return;
    
    for(auto it = _list.begin(); it != _list.end(); it++) {
        
        //if within range, add it
        if(it->first->distance(b->coordinate) < _rMax)
            it->second.push_back(b);
    }
}

void BoundaryBubbleNL::removeDynamicNeighbor(DynamicNeighbor* n) {
    
    //return if not a filament bead!
    Bubble* b;
    
    if(!(b = dynamic_cast<Bubble*>(n))) return;
    
    for(auto it = _list.begin(); it != _list.end(); it++) {
        
        auto bit = find(it->second.begin(), it->second.end(), b);
        if(bit != it->second.end()) it->second.erase(bit);
    }
}

void BoundaryBubbleNL::reset() {
    
    _list.clear();
    
    //loop through all neighbor keys
    for(auto boundary: BoundaryElement::getBoundaryElements())
        
        updateNeighbors(boundary);
}

vector<Bubble*> BoundaryBubbleNL::getNeighbors(BoundaryElement* be) {
    
    return _list[be];
}

//BUBBLE - BUBBLE

void BubbleBubbleNL::updateNeighbors(Bubble* bb) {
    
    //clear existing
    _list[bb].clear();
    
    //loop through beads, add as neighbor
    for (auto &bbo : Bubble::getBubbles()) {
        
        double dist = twoPointDistance(bb->coordinate, bbo->coordinate);
        
        if(bb->getID() <= bbo->getID()) continue;
        
        //If within range, add it
        if(dist < _rMax) _list[bb].push_back(bbo);
    }
}

void BubbleBubbleNL::addNeighbor(Neighbor* n) {
    
    //return if not a bubble!
    Bubble* bb;
    if(!(bb = dynamic_cast<Bubble*>(n))) return;
    
    //update neighbors
    updateNeighbors(bb);
}

void BubbleBubbleNL::removeNeighbor(Neighbor* n) {
    
    Bubble* bb;
    if(!(bb = dynamic_cast<Bubble*>(n))) return;
    
    _list.erase(bb);
    
    //remove from other lists
    for(auto it = _list.begin(); it != _list.end(); it++) {
        
        auto bit = find(it->second.begin(), it->second.end(), bb);
        if(bit != it->second.end()) it->second.erase(bit);
    }
}

void BubbleBubbleNL::reset() {
    
    _list.clear();
    
    //loop through all neighbor keys
    for(auto bb: Bubble::getBubbles())
        updateNeighbors(bb);
}

vector<Bubble*> BubbleBubbleNL::getNeighbors(Bubble* bb) {
    
    return _list[bb];
}

///BUBBLE - CYLINDER

void BubbleCylinderNL::updateNeighbors(Bubble* bb) {
    
    //clear existing
    _list[bb].clear();
    
    //loop through beads, add as neighbor
    for (auto &c : Cylinder::getCylinders()) {
        
        double dist = twoPointDistance(c->coordinate, bb->coordinate);
        
        //If within range, add it
        if(dist < _rMax) _list[bb].push_back(c);
    }
}

void BubbleCylinderNL::addNeighbor(Neighbor* n) {
    
    Bubble* bb; Cylinder* c;
    if((bb = dynamic_cast<Bubble*>(n))) {
        updateNeighbors(bb);
    }
    else if((c = dynamic_cast<Cylinder*>(n))) {
        
        for(auto it = _list.begin(); it != _list.end(); it++) {
            
            //if within range, add it
            if(twoPointDistance(it->first->coordinate, c->coordinate) < _rMax)
                it->second.push_back(c);
        }
    }
    else return;
}

void BubbleCylinderNL::removeNeighbor(Neighbor* n) {
    
    Bubble* bb; Cylinder* c;
    if((bb = dynamic_cast<Bubble*>(n))) {
        _list.erase(bb);
    }
    else if((c = dynamic_cast<Cylinder*>(n))) {
        for(auto it = _list.begin(); it != _list.end(); it++) {
            
            auto cit = find(it->second.begin(), it->second.end(), c);
            if(cit != it->second.end()) it->second.erase(cit);
        }
    }
    else return;
}

void BubbleCylinderNL::reset() {
    
    _list.clear();
    
    //loop through all neighbor keys
    for(auto bb: Bubble::getBubbles())
        updateNeighbors(bb);
}

vector<Cylinder*> BubbleCylinderNL::getNeighbors(Bubble* bb) {
    
    return _list[bb];
}

