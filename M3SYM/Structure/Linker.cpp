
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

#include "Linker.h"

#include "Bead.h"
#include "Cylinder.h"
#include "ChemRNode.h"

#include "GController.h"
#include "SystemParameters.h"
#include "MathFunctions.h"

using namespace mathfunc;

Linker::Linker(Cylinder* c1, Cylinder* c2, short linkerType,
               double position1, double position2, bool creation)
    : _c1(c1), _c2(c2),
      _position1(position1), _position2(position2), _linkerType(linkerType) {
                    
    //Add to linker db
    LinkerDB::instance()->addLinker(this);
    _linkerID = LinkerDB::instance()->getLinkerID();
                                            
    _birthTime = tau();
        
    //Find compartment
    auto m1 =
        midPointCoordinate(_c1->getFirstBead()->coordinate,
                           _c1->getSecondBead()->coordinate, _position1);
    auto m2 =
        midPointCoordinate(_c2->getFirstBead()->coordinate,
                           _c2->getSecondBead()->coordinate, _position2);
          
    coordinate = midPointCoordinate(m1, m2, 0.5);

    try {_compartment = GController::getCompartment(coordinate);}
    catch (exception& e) { cout << e.what(); exit(EXIT_FAILURE);}
        
#ifdef CHEMISTRY
    _cLinker = unique_ptr<CLinker>(new CLinker(_compartment));
    _cLinker->setLinker(this);
        
    //Find species on cylinder that should be marked. If initialization,
    //this should be done. But, if this is because of a reaction callback,
    //it will have already been done.
    int pos1 = int(position1 * SystemParameters::Geometry().cylinderIntSize);
    int pos2 = int(position2 * SystemParameters::Geometry().cylinderIntSize);
    
    SpeciesLinker* sl1 =
        _c1->getCCylinder()->getCMonomer(pos1)->speciesLinker(linkerType);
    SpeciesLinker* sl2 =
        _c2->getCCylinder()->getCMonomer(pos2)->speciesLinker(linkerType);
        
    if(!creation) {
        SpeciesBound* se1 =
            _c1->getCCylinder()->getCMonomer(pos1)->speciesBound(0);
        sl1->up();
        se1->down();
        
        SpeciesBound* se2 =
            _c2->getCCylinder()->getCMonomer(pos2)->speciesBound(0);
        sl2->up();
        se2->down();
    }
        
    //attach this linker to the species
    _cLinker->setFirstSpecies(sl1);
    _cLinker->setSecondSpecies(sl2);
#endif
    
#ifdef MECHANICS
    _mLinker = unique_ptr<MLinker>(
        new MLinker(linkerType, position1, position2,
            _c1->getFirstBead()->coordinate, _c1->getSecondBead()->coordinate,
            _c2->getFirstBead()->coordinate, _c2->getSecondBead()->coordinate));
    _mLinker->setLinker(this);
#endif
}

Linker::~Linker() noexcept {
    //Remove from linker db
    LinkerDB::instance()->removeLinker(this);
}

void Linker::updatePosition() {
    
    //check if were still in same compartment
    auto m1 =
        midPointCoordinate(_c1->getFirstBead()->coordinate,
                           _c1->getSecondBead()->coordinate, _position1);
    auto m2 =
        midPointCoordinate(_c2->getFirstBead()->coordinate,
                           _c2->getSecondBead()->coordinate, _position2);
    
    coordinate = midPointCoordinate(m1, m2, 0.5);
    
    _mLinker->setLength(twoPointDistance(m1, m2));
    
    Compartment* c;
    try {c = GController::getCompartment(coordinate);}
    catch (exception& e) { cout << e.what(); exit(EXIT_FAILURE);}
    
    if(c != _compartment) {
        
        _compartment = c;
#ifdef CHEMISTRY
        SpeciesBound* firstSpecies = _cLinker->getFirstSpecies();
        SpeciesBound* secondSpecies = _cLinker->getSecondSpecies();
        
        CLinker* clone = _cLinker->clone(c);
        setCLinker(clone);
        
        _cLinker->setFirstSpecies(firstSpecies);
        _cLinker->setSecondSpecies(secondSpecies);
#endif
    }
}

/// @note - This function updates unbinding rates based on the
/// following exponential form:
///
///                 k = k_0 * exp(f * a / kT)
///
/// The function uses the motor's stretching force at the current
/// state to change this rate.

void Linker::updateReactionRates() {

    //current force on linker
    double force = _mLinker->stretchForce;
    
    //characteristic length
    double a = SystemParameters::DynamicRates().LDULength[_linkerType];
    
    //get all walking reactions
    Species* s1 = _cLinker->getFirstSpecies();
    Species* s2 = _cLinker->getSecondSpecies();
    
    for(auto r : s1->getRSpecies().reactantReactions()) {
        if(r->getReactionType() == ReactionType::LINKERUNBINDING) {
            
            float newRate = r->getBareRate() * exp( force * a / kT);
            r->setRate(newRate);
            r->getRNode()->activateReaction();
        }
    }
    for(auto r : s2->getRSpecies().reactantReactions()) {
        if(r->getReactionType() == ReactionType::LINKERUNBINDING) {
            
            float newRate = r->getBareRate() * exp( force * a / kT);
            r->setRate(newRate);
            r->getRNode()->activateReaction();
        }
    }
}

