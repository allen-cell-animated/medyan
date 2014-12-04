
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

#include "BranchPoint.h"

#include "Bead.h"
#include "Cylinder.h"

#include "GController.h"
#include "SystemParameters.h"
#include "MathFunctions.h"

using namespace mathfunc;

BranchPoint::BranchPoint(Cylinder* c1, Cylinder* c2, short branchType, double position, bool creation)
                         : _c1(c1), _c2(c2), _branchType(branchType), _position(position){
                             
    //Add to branch point db
    BranchPointDB::instance()->addBranchPoint(this);
    _branchID = BranchPointDB::instance()->getBranchID();
    
    _birthTime = tau();
    
    //Find compartment
    coordinate = MidPointCoordinate(_c1->getFirstBead()->coordinate, _c1->getSecondBead()->coordinate, _position);
    try {_compartment = GController::getCompartment(coordinate);}
    catch (exception& e) { cout << e.what(); exit(EXIT_FAILURE);}
    
#ifdef CHEMISTRY
    _cBranchPoint = unique_ptr<CBranchPoint>(new CBranchPoint(_compartment));
    _cBranchPoint->setBranchPoint(this);
#endif
    
#ifdef MECHANICS
    _mBranchPoint = unique_ptr<MBranchPoint>(
      new MBranchPoint(SystemParameters::Mechanics().BrStretchingK[branchType],
                       SystemParameters::Mechanics().BrStretchingL[branchType],
                       SystemParameters::Mechanics().BrBendingK[branchType],
                       SystemParameters::Mechanics().BrBendingTheta[branchType],
                       SystemParameters::Mechanics().BrTwistingK[branchType],
                       SystemParameters::Mechanics().BrTwistingPhi[branchType]));
    _mBranchPoint->setBranchPoint(this);
#endif
}

BranchPoint::~BranchPoint() {
    //Remove from branch point db
    BranchPointDB::instance()->removeBranchPoint(this);
}


void BranchPoint::updatePosition() {
    
    //Find compartment
    coordinate = MidPointCoordinate(_c1->getFirstBead()->coordinate, _c1->getSecondBead()->coordinate, _position);

    Compartment* c;
    try {c = GController::getCompartment(coordinate);}
    catch (exception& e) { cout << e.what(); exit(EXIT_FAILURE);}
    
    if(c != _compartment) {
        _compartment = c;
#ifdef CHEMISTRY
        CBranchPoint* clone = _cBranchPoint->clone(c);
        setCBranchPoint(clone);
#endif
    }
}

void BranchPoint::updateReactionRates() {}


