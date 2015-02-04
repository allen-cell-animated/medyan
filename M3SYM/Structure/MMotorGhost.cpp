
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

#include "MMotorGhost.h"

#include "SystemParameters.h"
#include "MathFunctions.h"

using namespace mathfunc;

MMotorGhost::MMotorGhost(int motorType, double position1, double position2,
                 const vector<double>& coord11, const vector<double>& coord12,
                 const vector<double>& coord21, const vector<double>& coord22) {
    
    _kStretch = SystemParameters::Mechanics().MStretchingK[motorType];
    
    auto m1 = midPointCoordinate(coord11, coord12, position1);
    auto m2 = midPointCoordinate(coord21, coord22, position2);
    _eqLength = twoPointDistance(m1, m2);
}

