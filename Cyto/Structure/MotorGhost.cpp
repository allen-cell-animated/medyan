//
//  MotorGhost.cpp
//  CytoMech
//
//  Created by Konstantin Popov on 4/16/14.
//  Copyright (c) 2014 Konstantin Popov. All rights reserved.
//

#include "MotorGhost.h"

#include "Bead.h"
#include "Cylinder.h"

#include "SystemParameters.h"
#include "MathFunctions.h"

#include "GController.h"

using namespace mathfunc;

MotorGhost::MotorGhost(Cylinder* c1, Cylinder* c2, short motorType, int motorID, double position1, double position2, bool creation)
    : _c1(c1), _c2(c2), _motorType(motorType), _motorID(motorID), _position1(position1), _position2(position2) {
    
    _birthTime = tau();
    
    ///Find compartment
    auto m1 = MidPointCoordinate(_c1->getFirstBead()->coordinate, _c1->getSecondBead()->coordinate, _position1);
    auto m2 = MidPointCoordinate(_c2->getFirstBead()->coordinate, _c2->getSecondBead()->coordinate, _position2);
    coordinate = MidPointCoordinate(m1, m2, 0.5);
    
    try {_compartment = GController::getCompartment(coordinate);}
    catch (exception& e) { cout << e.what(); exit(EXIT_FAILURE);}
    
#ifdef CHEMISTRY
    _cMotorGhost = unique_ptr<CMotorGhost>(new CMotorGhost(_compartment));
    _cMotorGhost->setMotorGhost(this);
    
    
    ///Find species on cylinder that should be marked. If initialization, this should be done. But,
    ///if this is because of a reaction callback, it will have already been done.
    
    int pos1 = int(position1 * SystemParameters::Geometry().cylinderIntSize);
    int pos2 = int(position2 * SystemParameters::Geometry().cylinderIntSize);
    
    SpeciesMotor* sm1 = _c1->getCCylinder()->getCMonomer(pos1)->speciesMotor(_motorType);
    SpeciesMotor* sm2 = _c2->getCCylinder()->getCMonomer(pos2)->speciesMotor(_motorType);
    
    if(!creation) {
        SpeciesBound* se1 = _c1->getCCylinder()->getCMonomer(pos1)->speciesBound(0);
        if(sm1->getN() != 0) {
            sm1->getRSpecies().up();
            se1->getRSpecies().down();
        }
        SpeciesBound* se2 = _c2->getCCylinder()->getCMonomer(pos2)->speciesBound(0);
        if(sm2->getN() != 0) {
            sm2->getRSpecies().up();
            se2->getRSpecies().down();
        }
    }
    
    ///attach this linker to the species
    _cMotorGhost->setFirstSpecies(sm1);
    _cMotorGhost->setSecondSpecies(sm2);
    
#endif
    
#ifdef MECHANICS
    _mMotorGhost = unique_ptr<MMotorGhost>(
            new MMotorGhost(SystemParameters::Mechanics().MStretchingK[motorType], position1, position2,
                _c1->getFirstBead()->coordinate, _c1->getSecondBead()->coordinate,
                _c2->getFirstBead()->coordinate, _c2->getSecondBead()->coordinate));
    _mMotorGhost->setMotorGhost(this);
#endif
    
}

MotorGhost::~MotorGhost() {
    
#ifdef CHEMISTRY
    ///Find species on cylinder that should be unmarked. This should be done if deleting because
    ///of a reaction callback, but needs to be done if deleting for other reasons
    
    int pos1 = int(_position1 * SystemParameters::Geometry().cylinderIntSize);
    int pos2 = int(_position2 * SystemParameters::Geometry().cylinderIntSize);
    
    SpeciesMotor* sm1 = _c1->getCCylinder()->getCMonomer(pos1)->speciesMotor(_motorType);
    SpeciesBound* se1 = _c1->getCCylinder()->getCMonomer(pos1)->speciesBound(0);
    SpeciesMotor* sm2 = _c2->getCCylinder()->getCMonomer(pos2)->speciesMotor(_motorType);
    SpeciesBound* se2 = _c2->getCCylinder()->getCMonomer(pos2)->speciesBound(0);
    
    if(sm1->getN() != 0) {
        sm1->getRSpecies().down();
        se1->getRSpecies().up();
    }
    if(sm2->getN() != 0) {
        sm2->getRSpecies().down();
        se2->getRSpecies().up();
    }
    
#endif //CHEMISTRY
    
}


void MotorGhost::updatePosition() {
    
    ///check if were still in same compartment
    auto m1 = MidPointCoordinate(_c1->getFirstBead()->coordinate, _c1->getSecondBead()->coordinate, _position1);
    auto m2 = MidPointCoordinate(_c2->getFirstBead()->coordinate, _c2->getSecondBead()->coordinate, _position2);
    coordinate = MidPointCoordinate(m1, m2, 0.5);
    
    Compartment* c;
    try {c = GController::getCompartment(coordinate);}
    catch (exception& e) { cout << e.what(); exit(EXIT_FAILURE);}
    
    if(c != _compartment) {
        
        _compartment = c;
#ifdef CHEMISTRY
        CMotorGhost* clone = _cMotorGhost->clone(c);
        
        clone->setFirstSpecies(_cMotorGhost->getFirstSpecies());
        clone->setSecondSpecies(_cMotorGhost->getSecondSpecies());
        
        setCMotorGhost(clone);
#endif
    }
}