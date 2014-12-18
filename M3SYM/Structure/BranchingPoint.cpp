
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

#include "BranchingPoint.h"

#include "Bead.h"
#include "Cylinder.h"

#include "GController.h"
#include "SystemParameters.h"
#include "MathFunctions.h"

using namespace mathfunc;

BranchingPoint::BranchingPoint(Cylinder* c1, Cylinder* c2,
                               short branchType, double position,
                               bool creation)
    : _c1(c1), _c2(c2), _branchType(branchType), _position(position){
                             
    //Add to branch point db
    BranchingPointDB::instance()->addBranchingPoint(this);
    _branchID = BranchingPointDB::instance()->getBranchID();
    
    _birthTime = tau();
    
    //Find compartment
    coordinate = midPointCoordinate(_c1->getFirstBead()->coordinate,
                                    _c1->getSecondBead()->coordinate, _position);
    try {_compartment = GController::getCompartment(coordinate);}
    catch (exception& e) { cout << e.what(); exit(EXIT_FAILURE);}
    
#ifdef CHEMISTRY
    _cBranchingPoint = unique_ptr<CBranchingPoint>(new CBranchingPoint(_compartment));
    _cBranchingPoint->setBranchingPoint(this);
#endif
    
#ifdef MECHANICS
    _mBranchingPoint = unique_ptr<MBranchingPoint>(new MBranchingPoint(branchType));
    _mBranchingPoint->setBranchingPoint(this);
#endif
}

BranchingPoint::~BranchingPoint() {
    //Remove from branch point db
    BranchingPointDB::instance()->removeBranchingPoint(this);
}


void BranchingPoint::updatePosition() {
    
    //Find compartment
    coordinate = midPointCoordinate(_c1->getFirstBead()->coordinate,
                                    _c1->getSecondBead()->coordinate, _position);

    Compartment* c;
    try {c = GController::getCompartment(coordinate);}
    catch (exception& e) { cout << e.what(); exit(EXIT_FAILURE);}
    
    if(c != _compartment) {
        _compartment = c;
#ifdef CHEMISTRY
        CBranchingPoint* clone = _cBranchingPoint->clone(c);
        setCBranchingPoint(clone);
#endif
    }
}

void BranchingPoint::updateReactionRates() {}


