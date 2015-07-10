
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
#include "Filament.h"

#include "SysParams.h"
#include "GController.h"
#include "MathFunctions.h"

using namespace mathfunc;

Database<Bead*> Bead::_beads;

Bead::Bead (vector<double> v, Filament* f, int positionFilament)

    : Trackable(true, false, true, false),
      coordinate(v), force(3, 0), forceAux(3, 0),
      _pFilament(f), _positionFilament(positionFilament), _birthTime(tau()) {
    
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
    catch (exception& e) {
        
        //print exception
        cout << e.what() << endl;
        
        printInfo();
        
        //also print filament info
        _pFilament->printInfo();
        
        //exit
        exit(EXIT_FAILURE);
    }

    if(c != _compartment) {
        //remove from old compartment, add to new
        _compartment->removeBead(this);
        _compartment = c;
        _compartment->addBead(this);
    }
}

void Bead::printInfo() {
    
    cout << endl;
    
    cout << "Bead: ptr = " << this << endl;
    cout << "Coordinates = " << coordinate[0] << ", " << coordinate[1] << ", " << coordinate[2] << endl;
    cout << "Previous coordinates in minimization = " << coordinateP[0] << ", " << coordinateP[1] << ", " << coordinateP[2] << endl;
    cout << "Previous coordinates before minimization = " << coordinateB[0] << ", " << coordinateB[1] << ", " << coordinateB[2] << endl;
    cout << "Forces = " << force[0] << ", " << force[1] << ", " << force[2] << endl;
    cout << "Auxiliary forces = " << forceAux[0] << ", " << forceAux[1] << ", " << forceAux[2] << endl;

    cout << "Position on filament = " << _positionFilament << endl;
    cout << "Birth time = " << _birthTime << endl;
    
    cout << endl;
}

