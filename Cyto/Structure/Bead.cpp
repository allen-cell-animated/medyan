//
//  Bead.cpp
//  CytoMech
//
//  Created by Konstantin Popov on 4/15/14.
//  Copyright (c) 2014 Konstantin Popov. All rights reserved.
//

#include "Bead.h"
#include "BoundaryElementDB.h"
#include "SystemParameters.h"

Bead::Bead (std::vector<double> v): coordinate(v), force(3, 0), forceAux(3, 0)
{
    ///Find compartment, add this bead
    try {_compartment = GController::getCompartment(v);}
    catch (std::exception& e) {std:: cout << e.what(); exit(EXIT_FAILURE);}
    _compartment->addBead(this);
    
    ///Add to list of boundary elements
    for (auto &be : *BoundaryElementDB::Instance(BoundaryElementDBKey())) {
        ///If within cutoff, add bead to this boundary element interaction list
        if(be->distance(v) <= SystemParameters::Boundaries().boundaryCutoff) {
            _boundaryElements.push_back(be);
            be->addBead(this);
        }
    }
}

void Bead::updateBoundaryElements() {
    
    ///First, update this bead's list
    for(auto it = _boundaryElements.begin(); it != _boundaryElements.end(); it++) {
        auto be = (*it);
        if (be == nullptr || be->distance(coordinate) > SystemParameters::Boundaries().boundaryCutoff) {
            _boundaryElements.erase(it);
            be->removeBead(this);
        }
    }
    
    ///add any new interacting boundary elements
    for(auto &be : *BoundaryElementDB::Instance(BoundaryElementDBKey())) {
        if(be->distance(coordinate) <= SystemParameters::Boundaries().boundaryCutoff) {
            _boundaryElements.push_back(be);
            be->addBead(this);
        }
    }
    
}