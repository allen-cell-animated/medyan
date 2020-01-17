
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

#include "MMotorGhost.h"

#include "SysParams.h"
#include "MathFunctions.h"

using namespace mathfunc;

MMotorGhost::MMotorGhost(int motorType, int numBoundHeads, floatingpoint position1, floatingpoint position2,
                        const vector<floatingpoint>& coord11, const vector<floatingpoint>& coord12,
                        const vector<floatingpoint>& coord21, const vector<floatingpoint>& coord22) {
    
    if(!SysParams::Mechanics().MStretchingK.empty())
        _kStretch = SysParams::Mechanics().MStretchingK[motorType] * numBoundHeads;
    
    auto m1 = midPointCoordinate(coord11, coord12, position1);
    auto m2 = midPointCoordinate(coord21, coord22, position2);
    _eqLength = twoPointDistance(m1, m2);
}

void MMotorGhost::initializerestart(int motorType, floatingpoint eqLength,
		floatingpoint numBoundHeads){

    if(SysParams::RUNSTATE){
        LOG(ERROR) << "initializerestart Function from MMotorGhost class can only be "
                      "called "
                      "during restart phase. Exiting.";
        throw std::logic_error("Illegal function call pattern");
    }

    if(numBoundHeads > 0) {
        if(!SysParams::Mechanics().MStretchingK.empty())
            _kStretch = SysParams::Mechanics().MStretchingK[motorType] * numBoundHeads;
    }
    _eqLength = eqLength;}

void MMotorGhost::setStretchingConstant(int motorType, floatingpoint numBoundHeads) {
    
    if(!SysParams::Mechanics().MStretchingK.empty())
        _kStretch = SysParams::Mechanics().MStretchingK[motorType] * numBoundHeads;
}

