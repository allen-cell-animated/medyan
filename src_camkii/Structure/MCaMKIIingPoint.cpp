
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v3.1
//
//  Copyright (2015-2016)  Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the MEDYAN web page for more information:
//  http://www.medyan.org
//------------------------------------------------------------------

#include "MCaMKIIingPoint.h"

#include "SysParams.h"

MCaMKIIingPoint::MCaMKIIingPoint(int camkiiType) {
    
    //set parameters
    if(!SysParams::Mechanics().BrStretchingK.empty()) {
        _kStretch = SysParams::Mechanics().BrStretchingK[camkiiType];
        _eqLength = SysParams::Mechanics().BrStretchingL[camkiiType];
    }
    
    if(!SysParams::Mechanics().BrBendingK.empty()) {
        _kBend = SysParams::Mechanics().BrBendingK[camkiiType];
        _eqTheta = SysParams::Mechanics().BrBendingTheta[camkiiType];
    }
    
    if(!SysParams::Mechanics().BrDihedralK.empty())
        _kDihedr = SysParams::Mechanics().BrDihedralK[camkiiType];
 
    if(!SysParams::Mechanics().BrPositionK.empty())
        _kPosition = SysParams::Mechanics().BrPositionK[camkiiType];
}
