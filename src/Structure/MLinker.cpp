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
#include "Util/Io/Log.hpp"

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

void MLinker::initializerestart(floatingpoint eqLength){
    if(SysParams::RUNSTATE){
        LOG(ERROR) << "initializerestart Function from MLinker class can only be called "
                      "during restart phase. Exiting.";
        throw std::logic_error("Illegal function call pattern");
    }
    _eqLength = eqLength;}
