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

using namespace mathfunc;

MotorGhost::MotorGhost(Cylinder* pc1, Cylinder* pc2, short motorType, double position1, double position2){
    
    ///Find compartment
    auto m1 = MidPointCoordinate(_pc1->GetFirstBead()->coordinate, _pc1->GetSecondBead()->coordinate, _position1);
    auto m2 = MidPointCoordinate(_pc2->GetFirstBead()->coordinate, _pc2->GetSecondBead()->coordinate, _position2);
    auto position = MidPointCoordinate(m1, m2, 0.5);
    
    try {_compartment = GController::getCompartment(position);}
    catch (std::exception& e) {std:: cout << e.what(); exit(EXIT_FAILURE);}
    
#ifdef CHEMISTRY
    _cMotorGhost = std::unique_ptr<CMotorGhost>(new CMotorGhost(_compartment));
    _cMotorGhost->setMotorGhost(this);
    
    
    ///Find species on cylinder that should be marked. If initialization, this should be done. But,
    ///if this is because of a reaction callback, it will have already been done.
    
    int pos1 = int(position1 * SystemParameters::Geometry().cylinderSize);
    int pos2 = int(position2 * SystemParameters::Geometry().cylinderSize);
    
    SpeciesMotor* sm1 = _pc1->getCCylinder()->getCMonomer(pos1)->speciesMotor(_motorType);
    SpeciesBound* se1 = _pc1->getCCylinder()->getCMonomer(pos1)->speciesBound(0);
    SpeciesMotor* sm2 = _pc2->getCCylinder()->getCMonomer(pos2)->speciesMotor(_motorType);
    SpeciesBound* se2 = _pc2->getCCylinder()->getCMonomer(pos2)->speciesBound(0);
    
    if(sm1->getN() != 0) {
        sm1->getRSpecies().up();
        se1->getRSpecies().down();
    }
    if(sm2->getN() != 0) {
        sm2->getRSpecies().up();
        se2->getRSpecies().down();
    }
    
    ///attach this linker to the species
    sm1->setCBound(this->getCMotorGhost());
    sm2->setCBound(this->getCMotorGhost());
    
#endif
    
#ifdef MECHANICS
    _mMotorGhost = std::unique_ptr<MMotorGhost>(
            new MMotorGhost(SystemParameters::Mechanics().MStretchingK[motorType], position1, position2,
                pc1->GetFirstBead()->coordinate, pc1->GetSecondBead()->coordinate,
                pc2->GetFirstBead()->coordinate, pc2->GetSecondBead()->coordinate));
    _mMotorGhost->setMotorGhost(this);
#endif
    
}

MotorGhost::~MotorGhost() {
    
#ifdef CHEMISTRY
    ///Find species on cylinder that should be unmarked. This should be done if deleting because
    ///of a reaction callback, but needs to be done if deleting for other reasons
    
    int pos1 = int(_position1 * SystemParameters::Geometry().cylinderSize);
    int pos2 = int(_position2 * SystemParameters::Geometry().cylinderSize);
    
    SpeciesMotor* sm1 = _pc1->getCCylinder()->getCMonomer(pos1)->speciesMotor(_motorType);
    SpeciesBound* se1 = _pc1->getCCylinder()->getCMonomer(pos1)->speciesBound(0);
    SpeciesMotor* sm2 = _pc2->getCCylinder()->getCMonomer(pos2)->speciesMotor(_motorType);
    SpeciesBound* se2 = _pc2->getCCylinder()->getCMonomer(pos2)->speciesBound(0);
    
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
    auto m1 = MidPointCoordinate(_pc1->GetFirstBead()->coordinate, _pc1->GetSecondBead()->coordinate, _position1);
    auto m2 = MidPointCoordinate(_pc2->GetFirstBead()->coordinate, _pc2->GetSecondBead()->coordinate, _position2);
    auto position = MidPointCoordinate(m1, m2, 0.5);
    
    Compartment* c;
    try {c = GController::getCompartment(position);}
    catch (std::exception& e) {std:: cout << e.what(); exit(EXIT_FAILURE);}
    
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