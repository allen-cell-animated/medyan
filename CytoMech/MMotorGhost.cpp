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
   
    _kStretch = stretchConst;
    position1 = pos1;
    position2 = pos2;
    
     _eqLength = TwoPointDistance(MidPointCoordinate(_pc1->GetFirstBead()->coordinate, _pc1->GetSecondBead()->coordinate, position1 ), MidPointCoordinate(_pc2->GetFirstBead()->coordinate, _pc2->GetSecondBead()->coordinate, position2 ) );
    
    _pNetwork = pn;
    
}
//Public mechanical interfaces which call private "potentials" and force expressions to calculate bonded interactions:

double MotorGhost::CopmuteEnergy(){
    
    double Us = EnergyHarmonicStretching(_pc1->GetFirstBead(), _pc1->GetSecondBead(), _pc2->GetFirstBead(), _pc2->GetSecondBead());
    
    return Us;
    
}

double MotorGhost::CopmuteEnergy(double d){
    
    double Us = EnergyHarmonicStretching(_pc1->GetFirstBead(), _pc1->GetSecondBead(), _pc2->GetFirstBead(), _pc2->GetSecondBead(), d);
    
    return Us;
}



void MotorGhost::CopmuteForce(){
    
    ForceHarmonicStretching(_pc1->GetFirstBead(), _pc1->GetSecondBead(), _pc2->GetFirstBead(), _pc2->GetSecondBead());
}



void MotorGhost::CopmuteForceAux(){
    
    ForceHarmonicStretchingAux(_pc1->GetFirstBead(), _pc1->GetSecondBead(), _pc2->GetFirstBead(), _pc2->GetSecondBead()); //Stretching.
}
