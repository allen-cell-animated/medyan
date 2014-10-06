//
//  Linker.cpp
//  Cyto
//
//  Created by James Komianos on 10/6/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#include "Linker.h"
#include "Bead.h"
#include "SystemParameters.h"

Linker::Linker(Cylinder* pc1, Cylinder* pc2, Compartment* c, double position1, double position2) {
    
    _pc1 = pc1;
    _pc2 = pc2;
    
#ifdef CHEMISTRY
    _cLinker = std::unique_ptr<CLinker>(new CLinker(c));
#endif
    
#ifdef MECHANICS
    _mLinker = std::unique_ptr<MLinker>(new MLinker(SystemParameters::Mechanics().LStretchingK, position1, position2,
                                        pc1->GetFirstBead()->coordinate, pc1->GetSecondBead()->coordinate,
                                        pc2->GetFirstBead()->coordinate, pc2->GetSecondBead()->coordinate));
#endif
    

}

void Linker::updatePosition() {
#ifdef CHEMISTRY
    
    ///check if were still in same compartment
    auto m1 = MidPointCoordinate(_pc1->GetFirstBead()->coordinate, _pc1->GetSecondBead()->coordinate, _position1);
    auto m2 = MidPointCoordinate(_pc2->GetFirstBead()->coordinate, _pc2->GetSecondBead()->coordinate, _position2);
    auto position = MidPointCoordinate(m1, m2, 0.5);
    
    Compartment* c;
    try {c = GController::getCompartment(position);}
    catch (std::exception& e) {std:: cout << e.what(); exit(EXIT_FAILURE);}
    
    if(c != _cLinker->getCompartment()) {
        CLinker* clone = _cLinker->clone(c);
        setCLinker(clone);
    }
#endif
}
