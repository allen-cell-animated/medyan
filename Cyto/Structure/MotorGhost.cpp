//
//  MotorGhost.cpp
//  CytoMech
//
//  Created by Konstantin Popov on 4/16/14.
//  Copyright (c) 2014 Konstantin Popov. All rights reserved.
//

#include "MotorGhost.h"
#include "Bead.h"

using namespace mathfunc;

MotorGhost::MotorGhost(Cylinder* pc1, Cylinder* pc2, double stretchConst, double pos1, double pos2){
    
    _pc1 = pc1;
    _pc2 = pc2;
    _kStretch = SystemParameters::Mechanics().MStretchingK;
    _position1 = pos1;
    _position2 = pos2;
    
    auto m1 = MidPointCoordinate(_pc1->getMCylinder()->GetFirstBead()->coordinate, _pc1->getMCylinder()->GetSecondBead()->coordinate, _position1 );
    auto m2 = MidPointCoordinate(_pc2->getMCylinder()->GetFirstBead()->coordinate, _pc2->getMCylinder()->GetSecondBead()->coordinate, _position2 );
    
     _eqLength = TwoPointDistance(m1, m2);
}
