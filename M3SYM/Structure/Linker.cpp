
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

#include <cmath>

#include "Linker.h"

#include "Bead.h"
#include "Cylinder.h"
#include "ChemRNode.h"

#include "GController.h"
#include "SysParams.h"
#include "MathFunctions.h"

using namespace mathfunc;

vector<LinkerRateChanger*> Linker::_unbindingChangers;

Linker::Linker(Cylinder* c1, Cylinder* c2, short linkerType,
               double position1, double position2, bool creation)

    : _c1(c1), _c2(c2),
      _position1(position1), _position2(position2), _linkerType(linkerType) {
                    
    //Add to linker db
    LinkerDB::instance()->addLinker(this);
          
    _linkerID = LinkerDB::instance()->getLinkerID();
    _birthTime = tau();
        
    auto c1b1 = _c1->getFirstBead()->coordinate;
    auto c1b2 = _c1->getSecondBead()->coordinate;
    auto c2b1 = _c2->getFirstBead()->coordinate;
    auto c2b2 = _c2->getSecondBead()->coordinate;
          
    //Find compartment
    auto m1 = midPointCoordinate(c1b1, c1b2, _position1);
    auto m2 = midPointCoordinate(c2b1, c2b2, _position2);
          
    coordinate = midPointCoordinate(m1, m2, 0.5);

    try {_compartment = GController::getCompartment(coordinate);}
    catch (exception& e) { cout << e.what(); exit(EXIT_FAILURE);}
          
    int pos1 = int(position1 * SysParams::Geometry().cylinderIntSize);
    int pos2 = int(position2 * SysParams::Geometry().cylinderIntSize);
        
#ifdef CHEMISTRY
    _cLinker = unique_ptr<CLinker>(
        new CLinker(linkerType, _compartment,
                    _c1->getCCylinder(), _c2->getCCylinder(), pos1, pos2));
    _cLinker->setLinker(this);
        
#endif
    
#ifdef MECHANICS
    _mLinker = unique_ptr<MLinker>(
        new MLinker(linkerType, position1, position2,
                    c1b1, c1b2, c2b1, c2b2));
    _mLinker->setLinker(this);
#endif
}

Linker::~Linker() noexcept {
    
    //Remove from linker db
    LinkerDB::instance()->removeLinker(this);
    
}

void Linker::updatePosition() {
    
#ifdef CHEMISTRY
    //update ccylinders
    _cLinker->setFirstCCylinder(_c1->getCCylinder());
    _cLinker->setSecondCCylinder(_c2->getCCylinder());
    
#endif
    
    //check if were still in same compartment
    auto c1b1 = _c1->getFirstBead()->coordinate;
    auto c1b2 = _c1->getSecondBead()->coordinate;
    auto c2b1 = _c2->getFirstBead()->coordinate;
    auto c2b2 = _c2->getSecondBead()->coordinate;

    auto m1 = midPointCoordinate(c1b1, c1b2, _position1);
    auto m2 = midPointCoordinate(c2b1, c2b2, _position2);
    
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

/// @note - The function uses the motor's stretching force at
/// the current state to change this rate. Does not consider
/// compression forces, only stretching.

void Linker::updateReactionRates() {

    //if no rate changers were defined, skip
    if(_unbindingChangers.empty()) return;
    
    //current force on linker
    double force = max(0.0, _mLinker->stretchForce);
    
    //get the unbinding reaction
    ReactionBase* offRxn = _cLinker->getOffReaction();

    //change the rate
    float newRate =
        _unbindingChangers[_linkerType]->
        changeRate(offRxn->getBareRate(), force);
    
    offRxn->setRate(newRate);
    offRxn->activateReaction();
}

