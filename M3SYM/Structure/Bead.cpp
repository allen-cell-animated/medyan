
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
#include "Bead.h"

#include "Compartment.h"
#include "GController.h"

#include "MathFunctions.h"
#include "SystemParameters.h"

using namespace mathfunc;

Bead::Bead (vector<double> v, int positionFilament): _positionFilament(positionFilament),
                                                     coordinate(v), coordinateAux(v),
                                                     force(3, 0), forceAux(3, 0) {
    //set birth time
    _birthTime = tau();
    
    //Find compartment, add this bead
    try {_compartment = GController::getCompartment(v);}
    catch (exception& e) {cout << e.what(); exit(EXIT_FAILURE);}
    _compartment->addBead(this);
}

Bead::~Bead() {
    //remove from compartment
    _compartment->removeBead(this);
}

void Bead::updatePosition() {
    
    //Update the compartment
    Compartment* c;
    try {c = GController::getCompartment(coordinate);}
    catch (exception& e) {cout << e.what(); exit(EXIT_FAILURE);}
    
    if(c != _compartment) {
        //remove from old compartment, add to new
        _compartment->removeBead(this);
        _compartment = c;
        _compartment->addBead(this);
    }
}