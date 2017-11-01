
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

#include "MMotorGhost.h"

#include "SysParams.h"
#include "MathFunctions.h"

using namespace mathfunc;

MMotorGhost::MMotorGhost(int motorType, int numBoundHeads, double position1, double position2,
                        const vector<double>& coord11, const vector<double>& coord12,
                        const vector<double>& coord21, const vector<double>& coord22) {
    
    if(!SysParams::Mechanics().MStretchingK.empty())
        _kStretch = SysParams::Mechanics().MStretchingK[motorType] * numBoundHeads;
    
    auto m1 = midPointCoordinate(coord11, coord12, position1);
    auto m2 = midPointCoordinate(coord21, coord22, position2);
    _eqLength = twoPointDistance(m1, m2);
}

void MMotorGhost::setStretchingConstant(int motorType, double numBoundHeads) {
    
    //to account for bipolar nature of effective spring, divide by 2.
    if(!SysParams::Mechanics().MStretchingK.empty())
        _kStretch = SysParams::Mechanics().MStretchingK[motorType] * numBoundHeads / 2.0;
}

