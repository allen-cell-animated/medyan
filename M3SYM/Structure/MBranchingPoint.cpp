
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

#include "MBranchingPoint.h"

#include "SystemParameters.h"

MBranchingPoint::MBranchingPoint(int branchType) {
    
    //set parameters
    if(SystemParameters::Mechanics().BrStretchingK.size() != 0) {
        _kStretch = SystemParameters::Mechanics().BrStretchingK[branchType];
        _eqLength = SystemParameters::Mechanics().BrStretchingL[branchType];
    }
    if(SystemParameters::Mechanics().BrBendingK.size() != 0) {
        _kBend = SystemParameters::Mechanics().BrBendingK[branchType];
        _eqTheta = SystemParameters::Mechanics().BrBendingTheta[branchType];
    }
    if(SystemParameters::Mechanics().BrDihedralK.size() != 0)
        _kDihedr = SystemParameters::Mechanics().BrDihedralK[branchType];
}