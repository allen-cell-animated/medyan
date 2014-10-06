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

Bead::Bead (std::vector<double> v, int ID): coordinate(v), force(3, 0), forceAux(3, 0), _ID(ID)
{
    ///set birth time
    _birthTime = tau();
    
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

void Bead::updatePosition() {
    
    ///First, update this bead's list
    
    std::vector<BoundaryElement*> _beToRemove;
    for(auto &be : _boundaryElements) {
        if (be->distance(coordinate) > SystemParameters::Boundaries().boundaryCutoff) {
            _beToRemove.push_back(be); be->removeBead(this);
        }
    }
    for(auto &be : _beToRemove) 
        _boundaryElements.erase(std::find(_boundaryElements.begin(), _boundaryElements.end(), be));
    
    
    ///add any new interacting boundary elements
    for(auto &be : *BoundaryElementDB::Instance(BoundaryElementDBKey())) {
        if(be->distance(coordinate) <= SystemParameters::Boundaries().boundaryCutoff) {
            _boundaryElements.push_back(be);
            be->addBead(this);
        }
    }
    
}