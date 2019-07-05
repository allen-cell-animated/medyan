//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v4.0
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

MLinker::MLinker(int linkerType, floatingpoint position1, floatingpoint position2,
                 const vector<floatingpoint>& coord11, const vector<floatingpoint>& coord12,
                 const vector<floatingpoint>& coord21, const vector<floatingpoint>& coord22) {
    
    if(!SysParams::Mechanics().LStretchingK.empty())
        _kStretch = SysParams::Mechanics().LStretchingK[linkerType];
    
    auto m1 = midPointCoordinate(coord11, coord12, position1);
    auto m2 = midPointCoordinate(coord21, coord22, position2);
    _eqLength = twoPointDistance(m1, m2);
}
