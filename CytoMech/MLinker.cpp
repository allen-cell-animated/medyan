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



Linker::Linker(Network* pn ,Cylinder* pc1, Cylinder* pc2, double stretchConst, double position1, double position2){
    
    
    _pNetwork = pn;
    _pc1 = pc1;
    _pc2 = pc2;
    _kStretch = stretchConst;
    _position1 = position1;
    _position2 = position2;
    
    _eqLength = TwoPointDistance(MidPointCoordinate(_pc1->GetFirstBead()->coordinate, _pc1->GetSecondBead()->coordinate, position1 ), MidPointCoordinate(_pc2->GetFirstBead()->coordinate, _pc2->GetSecondBead()->coordinate, position2 ) );
    
    // For Tests:
    std::cout<<"Link Bead1 coord"<<_pc1->GetFirstBead()->coordinate[0]<<"  "<<_pc1->GetFirstBead()->coordinate[1]<<"  "<<_pc1->GetFirstBead()->coordinate[2]<<std::endl;
    std::cout<<"Link Bead2 coord"<<_pc2->GetFirstBead()->coordinate[0]<<"  "<<_pc2->GetFirstBead()->coordinate[1]<<"  "<<_pc2->GetFirstBead()->coordinate[2]<<std::endl;
    
}

double Linker::CopmuteEnergy(){
    
    double Us = EnergyHarmonicStretching(_pc1->GetFirstBead(), _pc2->GetFirstBead());
    
    return Us;
    
}

double Linker::CopmuteEnergy(double d){
    
    double Us = EnergyHarmonicStretching(_pc1->GetFirstBead(), _pc2->GetFirstBead(), d);
    
    return Us;
}



void Linker::CopmuteForce(){
    
    ForceHarmonicStretching(_pc1->GetFirstBead(), _pc2->GetFirstBead());
}



void Linker::CopmuteForceAux(){
    
    ForceHarmonicStretchingAux(_pc1->GetFirstBead(), _pc2->GetFirstBead()); //Stretching.
}
