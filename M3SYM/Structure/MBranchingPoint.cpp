
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
    if(!SysParams::Mechanics().BrStretchingK.empty()) {
        _kStretch = SysParams::Mechanics().BrStretchingK[branchType];
        _eqLength = SysParams::Mechanics().BrStretchingL[branchType];
    }
    
    if(!SysParams::Mechanics().BrBendingK.empty()) {
        _kBend = SysParams::Mechanics().BrBendingK[branchType];
        _eqTheta = SysParams::Mechanics().BrBendingTheta[branchType];
    }
    
    if(!SysParams::Mechanics().BrDihedralK.empty())
        _kDihedr = SysParams::Mechanics().BrDihedralK[branchType];
 
    if(!SysParams::Mechanics().BrPositionK.empty())
        _kPosition = SysParams::Mechanics().BrPositionK[branchType];
}