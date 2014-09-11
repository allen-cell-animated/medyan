//
//  Bead.cpp
//  CytoMech
//
//  Created by Konstantin Popov on 4/15/14.
//  Copyright (c) 2014 Konstantin Popov. All rights reserved.
//

#include "Bead.h"

Bead::Bead (std::vector<double> v): coordinate(v), force(3, 0), forceAux(3, 0)
{
    ///Find compartment, add this bead
    _compartment = GController::getCompartment(v);
    _compartment->addBead(this);
    
    ///Add to list of boundary elements in compartment
    for (auto &be : _compartment->getBoundaryElements()) {
        ///If within cutoff, add bead to this boundary element interaction list
        if(TwoPointDistance(be->coords(), coordinate) <= SystemParameters::Boundaries().interactionCutoff) {
            _boundaryElements.push_back(be);
            be->addBead(this);
        }
    }
}

void Bead::updateBoundaryElements() {
    
    ///First, update this bead's list
    for(auto it = _boundaryElements.begin(); it != _boundaryElements.end(); it++) {
        auto be = (*it);
        if (be == nullptr || TwoPointDistance(be->coords(), coordinate) > SystemParameters::Boundaries().interactionCutoff) {
            _boundaryElements.erase(it);
            be->removeBead(this);
        }
    }
    
    ///Check compartment, add any new interacting boundary elements
    for(auto &be : _compartment->getBoundaryElements()) {
        if(TwoPointDistance(be->coords(), coordinate) <= SystemParameters::Boundaries().interactionCutoff) {
            _boundaryElements.push_back(be);
            be->addBead(this);
        }
    }
    
}