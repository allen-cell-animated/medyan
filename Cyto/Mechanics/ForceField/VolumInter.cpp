//
//  VolumInter.cpp
//  CytoMech
//
//  Created by Konstantin Popov on 6/25/14.
//  Copyright (c) 2014 Konstantin Popov. All rights reserved.
//

#include "VolumInter.h"
#include "Bead.h"
#include "Cylinder.h"
//
//VolumInter::VolumInter(Network* pn, Cylinder* pc1, Cylinder* pc2, double repusionConst)
//{
//    
//    _pc1 = pc1;
//    _pc2 = pc2;
//    _pNetwork = pn;
//    _kRepuls = repusionConst;
//}
//
//
////Public mechanical interfaces which call private "potentials" and force expressions to calculate bonded interactions:
//
//double VolumInter::ComputeEnergy(){
//    
//    ///Iterate over some neighbour list and compute volume interactions.
//    
//    double U = 0;
//    
//    
//    U += EnergyRepulsion(_pc1->GetFirstBead(), _pc1->GetSecondBead(), _pc2->GetFirstBead(), _pc2->GetSecondBead());
//    
//    return U;
//    
//}
//
//double VolumInter::ComputeEnergy(double d){
//    
//    ///Iterate over some neighbour list and compute volume interactions.
//    
//    double U = 0;
//    
//    
//    U += EnergyRepulsion(_pc1->GetFirstBead(), _pc1->GetSecondBead(), _pc2->GetFirstBead(), _pc2->GetSecondBead(), d);
//    
//    return U;
//    
//}
//
///// Interface that iterates over neighbourlist of cylinders and compute pair repulsive forces ( on 4 beads);
//void VolumInter:: CopmuteForce(){
//    
//    ///Iterate over some neighbour list and compute volume interactions.
//    
//    
//        
//   // ForceRepulsion(_pc1->GetFirstBead(), _pc1->GetSecondBead(), _pc2->GetFirstBead(), _pc2->GetSecondBead() );
//    
//    
//        
//}
//
//void VolumInter:: CopmuteForceAux(){
//    
//    ///Iterate over some neighbour list and compute volume interactions.
//        
//  //  ForceRepulsionAux(_pc1->GetFirstBead(), _pc1->GetSecondBead(), _pc2->GetFirstBead(), _pc2->GetSecondBead() );
//
//}
//
//
