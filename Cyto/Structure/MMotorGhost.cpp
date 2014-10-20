//
//  MMotorGhost.cpp
//  Cyto
//
//  Created by James Komianos on 10/20/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#include "MMotorGhost.h"
#include "SystemParameters.h"
#include "MathFunctions.h"

using namespace mathfunc;

MMotorGhost::MMotorGhost(double stretchConst, double position1, double position2,
                 const std::vector<double>& coord11, const std::vector<double>& coord12,
                 const std::vector<double>& coord21, const std::vector<double>& coord22) {
    
    _kStretch = stretchConst;
    
    auto m1 = MidPointCoordinate(coord11, coord12, position1);
    auto m2 = MidPointCoordinate(coord21, coord22, position2);
    _eqLength = TwoPointDistance(m1, m2);
}

