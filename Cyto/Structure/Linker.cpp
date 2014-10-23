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

Linker::Linker(Cylinder* pc1, Cylinder* pc2, short linkerType, double position1, double position2, bool creation) :
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
        
        
    ///Find species on cylinder that should be marked. If initialization, this should be done. But,
    ///if this is because of a reaction callback, it will have already been done.
        
    int pos1 = int(position1 * SystemParameters::Geometry().cylinderIntSize);
    int pos2 = int(position2 * SystemParameters::Geometry().cylinderIntSize);
    
    SpeciesLinker* sl1 = pc1->getCCylinder()->getCMonomer(pos1)->speciesLinker(linkerType);
    SpeciesLinker* sl2 = pc2->getCCylinder()->getCMonomer(pos2)->speciesLinker(linkerType);
        
    if(!creation) {
        if(sl1->getN() != 1) {
            SpeciesBound* se1 = _pc1->getCCylinder()->getCMonomer(pos1)->speciesBound(0);
            sl1->getRSpecies().up();
            se1->getRSpecies().down();
        }
        if(sl2->getN() != 1) {
            SpeciesBound* se2 = _pc2->getCCylinder()->getCMonomer(pos2)->speciesBound(0);
            sl2->getRSpecies().up();
            se2->getRSpecies().down();
        }
    }
        
    ///attach this linker to the species
    _cLinker->setFirstSpecies(sl1);
    _cLinker->setSecondSpecies(sl2);
#endif
    
#ifdef MECHANICS
    _mLinker = std::unique_ptr<MLinker>(
                    new MLinker(SystemParameters::Mechanics().LStretchingK[linkerType], position1, position2,
                                        pc1->GetFirstBead()->coordinate, pc1->GetSecondBead()->coordinate,
                                        pc2->GetFirstBead()->coordinate, pc2->GetSecondBead()->coordinate));
    _mLinker->setLinker(this);
#endif
}

Linker::~Linker() {
    
#ifdef CHEMISTRY
    ///Find species on cylinder that should be unmarked. This should be done if deleting because
    ///of a reaction callback, but needs to be done if deleting for other reasons

    int pos1 = int(_position1 * SystemParameters::Geometry().cylinderIntSize);
    int pos2 = int(_position2 * SystemParameters::Geometry().cylinderIntSize);
    
    SpeciesLinker* sl1 = _pc1->getCCylinder()->getCMonomer(pos1)->speciesLinker(_linkerType);
    SpeciesBound* se1 = _pc1->getCCylinder()->getCMonomer(pos1)->speciesBound(0);
    SpeciesLinker* sl2 = _pc2->getCCylinder()->getCMonomer(pos2)->speciesLinker(_linkerType);
    SpeciesBound* se2 = _pc2->getCCylinder()->getCMonomer(pos2)->speciesBound(0);
    
    if(sl1->getN() != 0) {
        sl1->getRSpecies().down();
        se1->getRSpecies().up();
    }
    if(sl2->getN() != 0) {
        sl2->getRSpecies().down();
        se2->getRSpecies().up();
    }
    
#endif //CHEMISTRY
    
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
        
        clone->setFirstSpecies(_cLinker->getFirstSpecies());
        clone->setSecondSpecies(_cLinker->getSecondSpecies());
        
        setCLinker(clone);
#endif
    }
}
