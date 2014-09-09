//
//  MMotorGhost.cpp
//  CytoMech
//
//  Created by Konstantin Popov on 4/16/14.
//  Copyright (c) 2014 Konstantin Popov. All rights reserved.
//

#include "MMotorGhost.h"
#include "MBead.h"
#include "MCylinder.h"


using namespace mathfunc;


MotorGhost::MotorGhost(Network* pn, Cylinder* pc1, Cylinder* pc2, double stretchConst, double pos1, double pos2){
    
    _pc1 = pc1;
    _pc2 = pc2;
    _eqLength = SystemParameters::Mechanics().MStretchingL;
    _kStretch = SystemParameters::Mechanics().MStretchingK;
    _position1 = pos1;
    _position2 = pos2;
    
     _eqLength = TwoPointDistance(
        MidPointCoordinate(_pc1->GetFirstBead()->coordinate, _pc1->GetSecondBead()->coordinate, position1 ),
        MidPointCoordinate(_pc2->GetFirstBead()->coordinate, _pc2->GetSecondBead()->coordinate, position2 ));
}
