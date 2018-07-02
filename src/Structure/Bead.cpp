
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

#include "Bead.h"

#include "Compartment.h"
#include "Composite.h"

#include "SysParams.h"
#include "GController.h"
#include "MathFunctions.h"

using namespace mathfunc;

Database<Bead*> Bead::_beads;
Database<Bead*> Bead::_pinnedBeads;

Bead::Bead (vector<double> v, Composite* parent, int position)

    : Trackable(true),
      coordinate(v), coordinateP(v), coordinateB(v),
      force(3, 0), forceAux(3, 0), forceAuxP(3, 0),
      _position(position), _birthTime(tau()) {
    
    parent->addChild(unique_ptr<Component>(this));
          
    loadForcesP = vector<double>(SysParams::Geometry().cylinderNumMon[getType()], 0);
    loadForcesM = vector<double>(SysParams::Geometry().cylinderNumMon[getType()], 0);
    
    //Find compartment
    try {_compartment = GController::getCompartment(v);}
    catch (exception& e) {
        
        cout << e.what() << endl;
        
        printSelf();
        
        //also print parent info
        getParent()->printSelf();
        
        exit(EXIT_FAILURE);
    }
          
}

Bead::Bead(Composite* parent, int position)

    : Trackable(true),
    coordinate(3, 0), coordinateP(3, 0), coordinateB(3, 0),
    force(3, 0), forceAux(3, 0), forceAuxP(3, 0),  _position(position) {
    
    parent->addChild(unique_ptr<Component>(this));
}

void Bead::updatePosition() {
    
    try {GController::getCompartment(coordinate);}
    catch (exception& e) {
        
        //print exception
        cout << e.what() << endl;
        
        printSelf();
        
        //also print parent info
        getParent()->printSelf();
        
        //exit
        exit(EXIT_FAILURE);
    }
}

void Bead::printSelf() {
    
    cout << endl;
    
    cout << "Bead: ptr = " << this << endl;
    cout << "Coordinates = " << coordinate[0] << ", " << coordinate[1] << ", " << coordinate[2] << endl;
    cout << "Previous coordinates in minimization = " << coordinateP[0] << ", " << coordinateP[1] << ", " << coordinateP[2] << endl;
    cout << "Previous coordinates before minimization = " << coordinateB[0] << ", " << coordinateB[1] << ", " << coordinateB[2] << endl;
    cout << "Forces = " << force[0] << ", " << force[1] << ", " << force[2] << endl;
    cout << "Auxiliary forces = " << forceAux[0] << ", " << forceAux[1] << ", " << forceAux[2] << endl;

    cout << "Position on structure = " << _position << endl;
    cout << "Birth time = " << _birthTime << endl;
    
    cout << endl;
}

double Bead::getLoadForcesP() {
    
    if (lfip < 0)
        return loadForcesP[0];
        
    if (lfip >= SysParams::Geometry().cylinderNumMon[getType()])
        return loadForcesP[SysParams::Geometry().cylinderNumMon[getType()] - 1];
    
    else return loadForcesP[lfip];
}

double Bead::getLoadForcesM() {
    
    if (lfim < 0)
        return loadForcesM[0];
    
    if (lfim >= SysParams::Geometry().cylinderNumMon[getType()])
        return loadForcesM[SysParams::Geometry().cylinderNumMon[getType()] - 1];
    
    else return loadForcesM[lfim];
}

