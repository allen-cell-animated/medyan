//
//  Linker.cpp
//  Cyto
//
//  Created by James Komianos on 10/6/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#include "Linker.h"
#include "Bead.h"
#include "Cylinder.h"
#include "SystemParameters.h"

Linker::Linker(Cylinder* pc1, Cylinder* pc2, short linkerType, double position1, double position2) :
    _pc1(pc1), _pc2(pc2), _linkerType(linkerType), _position1(position1), _position2(position2) {

    ///Find compartment
    auto m1 = MidPointCoordinate(_pc1->GetFirstBead()->coordinate, _pc1->GetSecondBead()->coordinate, _position1);
    auto m2 = MidPointCoordinate(_pc2->GetFirstBead()->coordinate, _pc2->GetSecondBead()->coordinate, _position2);
    auto position = MidPointCoordinate(m1, m2, 0.5);

    try {_compartment = GController::getCompartment(position);}
    catch (std::exception& e) {std:: cout << e.what(); exit(EXIT_FAILURE);}
        
#ifdef CHEMISTRY
    _cLinker = std::unique_ptr<CLinker>(new CLinker(_compartment));
    _cLinker->setLinker(this);
#endif
    
#ifdef MECHANICS
    _mLinker = std::unique_ptr<MLinker>(
                    new MLinker(SystemParameters::Mechanics().LStretchingK[linkerType], position1, position2,
                                        pc1->GetFirstBead()->coordinate, pc1->GetSecondBead()->coordinate,
                                        pc2->GetFirstBead()->coordinate, pc2->GetSecondBead()->coordinate));
    _mLinker->setLinker(this);
#endif
    

}

void Linker::updatePosition() {
    
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
        CLinker* clone = _cLinker->clone(c);
        setCLinker(clone);
#endif
    }
}
