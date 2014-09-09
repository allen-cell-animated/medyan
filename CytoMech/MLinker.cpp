//
//  MLinker.cpp
//  CytoMech
//
//  Created by Konstantin Popov on 4/15/14.
//  Copyright (c) 2014 Konstantin Popov. All rights reserved.
//

#include "MLinker.h"
#include "MBead.h"
#include "MCylinder.h"

using namespace mathfunc;

Linker::Linker(Cylinder* pc1, Cylinder* pc2, double stretchConst, double position1, double position2) {
    
    _pc1 = pc1;
    _pc2 = pc2;
    _eqLength = SystemParameters::Mechanics().LStretchingL;
    _kStretch = SystemParameters::Mechanics().LStretchingK;
    _position1 = position1;
    _position2 = position2;
    
    _eqLength = TwoPointDistance(
        MidPointCoordinate(_pc1->GetFirstBead()->coordinate, _pc1->GetSecondBead()->coordinate, position1 ),
        MidPointCoordinate(_pc2->GetFirstBead()->coordinate, _pc2->GetSecondBead()->coordinate, position2 ));
}
