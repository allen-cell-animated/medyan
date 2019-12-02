
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
    if(!SysParams::Mechanics().CaMKIIStretchingK.empty()) {
        _kStretch = SysParams::Mechanics().CaMKIIStretchingK[camkiiType];
        _eqLength = SysParams::Mechanics().CaMKIIStretchingL[camkiiType];
    }
    
    if(!SysParams::Mechanics().CaMKIIBendingK.empty()) {
        _kBend = SysParams::Mechanics().CaMKIIBendingK[camkiiType];
        _eqTheta = SysParams::Mechanics().CaMKIIBendingTheta[camkiiType];
    }
    
    if(!SysParams::Mechanics().CaMKIIDihedralK.empty())
        _kDihedr = SysParams::Mechanics().CaMKIIDihedralK[camkiiType];
 
    if(!SysParams::Mechanics().CaMKIIPositionK.empty())
        _kPosition = SysParams::Mechanics().CaMKIIPositionK[camkiiType];
}
