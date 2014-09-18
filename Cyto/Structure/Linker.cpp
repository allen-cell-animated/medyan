//
//  Linker.cpp
//  CytoMech
//
//  Created by Konstantin Popov on 4/15/14.
//  Copyright (c) 2014 Konstantin Popov. All rights reserved.
//

#include "Linker.h"
#include "Bead.h"
#include "Cylinder.h"
#include "MathFunctions.h"
#include "SystemParameters.h"

using namespace mathfunc;

Linker::Linker(Cylinder* pc1, Cylinder* pc2, double stretchConst, double position1, double position2) {
    
    _pc1 = pc1;
    _pc2 = pc2;
    _kStretch = SystemParameters::Mechanics().LStretchingK;
    _position1 = position1;
    _position2 = position2;
    
    auto m1 = MidPointCoordinate(_pc1->getMCylinder()->GetFirstBead()->coordinate, _pc1->getMCylinder()->GetSecondBead()->coordinate, position1);
    auto m2 = MidPointCoordinate(_pc2->getMCylinder()->GetFirstBead()->coordinate, _pc2->getMCylinder()->GetSecondBead()->coordinate, position2);
    _eqLength = TwoPointDistance(m1, m2);
}
