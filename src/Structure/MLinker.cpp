
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v3.2
//
//  Copyright (2015-2018)  Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the MEDYAN web page for more information:
//  http://www.medyan.org
//------------------------------------------------------------------

#include "MLinker.h"

#include "SysParams.h"
#include "MathFunctions.h"

using namespace mathfunc;

MLinker::MLinker(int linkerType, double position1, double position2,
                 const vector<double>& coord11, const vector<double>& coord12,
                 const vector<double>& coord21, const vector<double>& coord22) {
    
    if(!SysParams::Mechanics().LStretchingK.empty())
        _kStretch = SysParams::Mechanics().LStretchingK[linkerType];
    
    auto m1 = midPointCoordinate(coord11, coord12, position1);
    auto m2 = midPointCoordinate(coord21, coord22, position2);
    _eqLength = twoPointDistance(m1, m2);
}
