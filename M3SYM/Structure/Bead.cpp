
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
#include "Bead.h"

#include "Compartment.h"
#include "NeighborListDB.h"

#include "GController.h"
#include "MathFunctions.h"
#include "SysParams.h"

using namespace mathfunc;

Bead::Bead (vector<double> v, int positionFilament)

    : coordinate(v), coordinateAux(v), force(3, 0), forceAux(3, 0),
      _positionFilament(positionFilament) {
          
    //add to bead db
    BeadDB::instance()->addBead(this);
                                                         
    //set birth time
    _birthTime = tau();
    
    //Find compartment, add this bead
    try {_compartment = GController::getCompartment(v);}
    catch (exception& e) {cout << e.what(); exit(EXIT_FAILURE);}
    _compartment->addBead(this);
                    
    //add to neighbor lists
    NeighborListDB::instance()->addDynamicNeighbor(this);
}

Bead::~Bead() {
    //remove from bead db
    BeadDB::instance()->removeBead(this);
    
    //remove from compartment
    _compartment->removeBead(this);
    
    //remove from neighbor lists
    NeighborListDB::instance()->removeDynamicNeighbor(this);
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