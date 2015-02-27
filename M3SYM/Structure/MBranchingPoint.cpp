
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

#include "SysParams.h"

MBranchingPoint::MBranchingPoint(int branchType) {
    
    //set parameters
    _kStretch = SysParams::Mechanics().BrStretchingK[branchType];
    _eqLength = SysParams::Mechanics().BrStretchingL[branchType];

    _kBend = SysParams::Mechanics().BrBendingK[branchType];
    _eqTheta = SysParams::Mechanics().BrBendingTheta[branchType];
    
    _kDihedr = SysParams::Mechanics().BrDihedralK[branchType];
 
    _kPosition = SysParams::Mechanics().BrPositionK[branchType];
}