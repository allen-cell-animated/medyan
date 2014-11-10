//
//  MLinker.cpp
//  CytoMech
//
//  Created by Konstantin Popov on 4/15/14.
//  Copyright (c) 2014 Konstantin Popov. All rights reserved.
//

#include "MLinker.h"
#include "SystemParameters.h"
#include "MathFunctions.h"

using namespace mathfunc;

MLinker::MLinker(double stretchConst, double position1, double position2,
                 const vector<double>& coord11, const vector<double>& coord12,
                 const vector<double>& coord21, const vector<double>& coord22) {
    
    _kStretch = stretchConst;
    
    auto m1 = MidPointCoordinate(coord11, coord12, position1);
    auto m2 = MidPointCoordinate(coord21, coord22, position2);
    _eqLength = TwoPointDistance(m1, m2);
}
