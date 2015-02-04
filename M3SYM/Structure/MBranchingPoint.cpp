
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

#include "MBranchingPoint.h"

#include "SystemParameters.h"

MBranchingPoint::MBranchingPoint(int branchType) {
    
    //set parameters
    _kStretch = SystemParameters::Mechanics().BrStretchingK[branchType];
    _eqLength = SystemParameters::Mechanics().BrStretchingL[branchType];

    _kBend = SystemParameters::Mechanics().BrBendingK[branchType];
    _eqTheta = SystemParameters::Mechanics().BrBendingTheta[branchType];
    
    _kDihedr = SystemParameters::Mechanics().BrDihedralK[branchType];
 
    _kPosition = SystemParameters::Mechanics().BrPositionK[branchType];
}