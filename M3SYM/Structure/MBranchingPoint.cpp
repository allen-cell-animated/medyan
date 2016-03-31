
//------------------------------------------------------------------
//  **M3SYM** - Simulation Package for the Mechanochemical
//              Dynamics of Active Networks, 3rd Generation
//
//  Copyright (2015)  Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the MEDYAN web page for more information:
//  http://www.medyan.org
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