
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

#include "MotorGhost.h"

#include "Bead.h"
#include "Cylinder.h"

#include "GController.h"
#include "SystemParameters.h"
#include "MathFunctions.h"

using namespace mathfunc;

MotorGhost::MotorGhost(Cylinder* c1, Cylinder* c2, short motorType, double position1, double position2, bool creation)
                               : _c1(c1), _c2(c2), _motorType(motorType), _position1(position1), _position2(position2) {
    
    //add to motor ghost db
    MotorGhostDB::instance()->addMotorGhost(this);
    _motorID = MotorGhostDB::instance()->getMotorID();
    _birthTime = tau();
    
    //Find compartment
    auto m1 = midPointCoordinate(_c1->getFirstBead()->coordinate, _c1->getSecondBead()->coordinate, _position1);
    auto m2 = midPointCoordinate(_c2->getFirstBead()->coordinate, _c2->getSecondBead()->coordinate, _position2);
    coordinate = midPointCoordinate(m1, m2, 0.5);
    
    try {_compartment = GController::getCompartment(coordinate);}
    catch (exception& e) { cout << e.what(); exit(EXIT_FAILURE);}
    
#ifdef CHEMISTRY
    _cMotorGhost = unique_ptr<CMotorGhost>(new CMotorGhost(_compartment));
    _cMotorGhost->setMotorGhost(this);
    
    //Find species on cylinder that should be marked. If initialization, this should be done. But,
    //if this is because of a reaction callback, it will have already been done.
    
    int pos1 = int(position1 * SystemParameters::Geometry().cylinderIntSize);
    int pos2 = int(position2 * SystemParameters::Geometry().cylinderIntSize);
    
    SpeciesMotor* sm1 = _c1->getCCylinder()->getCMonomer(pos1)->speciesMotor(_motorType);
    SpeciesMotor* sm2 = _c2->getCCylinder()->getCMonomer(pos2)->speciesMotor(_motorType);
    
    if(!creation) {
        SpeciesBound* se1 = _c1->getCCylinder()->getCMonomer(pos1)->speciesBound(0);
        sm1->getRSpecies().up();
        se1->getRSpecies().down();
        
        SpeciesBound* se2 = _c2->getCCylinder()->getCMonomer(pos2)->speciesBound(0);
        sm2->getRSpecies().up();
        se2->getRSpecies().down();
    }
    
    //attach this motor to the species
    _cMotorGhost->setFirstSpecies(sm1);
    _cMotorGhost->setSecondSpecies(sm2);
    
#endif
    
#ifdef MECHANICS
    _mMotorGhost = unique_ptr<MMotorGhost>(new MMotorGhost(motorType, position1, position2,
                                           _c1->getFirstBead()->coordinate, _c1->getSecondBead()->coordinate,
                                           _c2->getFirstBead()->coordinate, _c2->getSecondBead()->coordinate));
    _mMotorGhost->setMotorGhost(this);
#endif
    
}

MotorGhost::~MotorGhost() {
    
    //remove from motor ghost db
    MotorGhostDB::instance()->removeMotorGhost(this);
}

void MotorGhost::updatePosition() {
    
    //check if were still in same compartment
    auto m1 = midPointCoordinate(_c1->getFirstBead()->coordinate, _c1->getSecondBead()->coordinate, _position1);
    auto m2 = midPointCoordinate(_c2->getFirstBead()->coordinate, _c2->getSecondBead()->coordinate, _position2);
    coordinate = midPointCoordinate(m1, m2, 0.5);
    
    Compartment* c;
    try {c = GController::getCompartment(coordinate);}
    catch (exception& e) { cout << e.what(); exit(EXIT_FAILURE);}
    
    if(c != _compartment) {
        
        _compartment = c;
#ifdef CHEMISTRY
        SpeciesBound* firstSpecies = _cMotorGhost->getFirstSpecies();
        SpeciesBound* secondSpecies = _cMotorGhost->getSecondSpecies();
        
        CMotorGhost* clone = _cMotorGhost->clone(c);
        setCMotorGhost(clone);
        
        _cMotorGhost->setFirstSpecies(firstSpecies);
        _cMotorGhost->setSecondSpecies(secondSpecies);
#endif
    }
}

void MotorGhost::updateReactionRates() {}
