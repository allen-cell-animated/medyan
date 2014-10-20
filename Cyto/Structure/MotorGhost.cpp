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
#endif
    
#ifdef MECHANICS
    _mMotorGhost = std::unique_ptr<MMotorGhost>(
            new MMotorGhost(SystemParameters::Mechanics().MStretchingK[motorType], position1, position2,
                pc1->GetFirstBead()->coordinate, pc1->GetSecondBead()->coordinate,
                pc2->GetFirstBead()->coordinate, pc2->GetSecondBead()->coordinate));
    _mMotorGhost->setMotorGhost(this);
#endif
    
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
        setCMotorGhost(clone);
#endif
    }
}