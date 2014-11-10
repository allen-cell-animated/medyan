//
//  Bead.cpp
//  CytoMech
//
//  Created by Konstantin Popov on 4/15/14.
//  Copyright (c) 2014 Konstantin Popov. All rights reserved.
//

#include "Bead.h"

#include "BoundaryElementDB.h"
#include "MathFunctions.h"
#include "SystemParameters.h"

using namespace mathfunc;

Bead::Bead (vector<double> v, int ID): coordinate(v), coordinateAux(v), force(3, 0), forceAux(3, 0), _ID(ID)
{
    ///set birth time
    _birthTime = tau();
    
    ///Find compartment, add this bead
    try {_compartment = GController::getCompartment(v);}
    catch (exception& e) {cout << e.what(); exit(EXIT_FAILURE);}
    _compartment->addBead(this);
    
    ///Add to list of boundary elements
    for (auto &be : *BoundaryElementDB::instance(BoundaryElementDBKey())) {
        ///If within cutoff, add bead to this boundary element interaction list
        if(be->distance(v) <= SystemParameters::Boundaries().boundaryCutoff) {
            addBoundaryElement(be);
            be->addBead(this);
        }
    }
}

Bead::~Bead() {
    ///remove from compartment
    _compartment->removeBead(this);
    
    ///remove from boundary elements
    for(auto &be : _boundaryElements) be->removeBead(this);
}

void Bead::updatePosition() {
    
    ///Update the compartment
    
    Compartment* c;
    try {c = GController::getCompartment(coordinate);}
    catch (exception& e) {cout << e.what(); exit(EXIT_FAILURE);}
    
    if(c != _compartment) {
        ///remove from old compartment, add to new
        _compartment->removeBead(this);
        _compartment = c;
        _compartment->addBead(this);
    }
    
    ///Update this bead's list
    vector<BoundaryElement*> _beToRemove;
    for(auto &be : _boundaryElements) {
        if (be->distance(coordinate) > SystemParameters::Boundaries().boundaryCutoff) {
            _beToRemove.push_back(be); be->removeBead(this);
        }
    }
    for(auto &be : _beToRemove)  removeBoundaryElement(be);

    ///add any new interacting boundary elements
    for(auto &be : *BoundaryElementDB::instance(BoundaryElementDBKey())) {
        if(be->distance(coordinate) <= SystemParameters::Boundaries().boundaryCutoff) {
            addBoundaryElement(be);
            be->addBead(this);
        }
    }
    
}