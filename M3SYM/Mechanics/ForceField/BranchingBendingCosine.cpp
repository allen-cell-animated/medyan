//
//  BranchingBendingCosine.cpp
//  M3SYM
//
//  Created by Konstantin Popov on 12/3/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#include "BranchingBendingCosine.h"

#include "MathFunctions.h"
#include "Bead.h"

#include <math.h>

using namespace mathfunc;

double BranchingBendingCosine::energy(Bead* b1, Bead* b2, Bead* b3, Bead* b4, double k_bend, double phi0){
    
    
    
    double L1 = sqrt(ScalarProduct(b1->coordinate, b2->coordinate, b1->coordinate, b2->coordinate));
    double L2 = sqrt(ScalarProduct(b3->coordinate, b4->coordinate, b3->coordinate, b4->coordinate));
    
    double L1L2 = L1*L2;
    double l1l2 = ScalarProduct(b1->coordinate, b2->coordinate, b3->coordinate, b4->coordinate);
    
    double phi = acos(l1l2 / L1L2);
    double dPhi = phi-phi0;
    
    return k_bend * ( 1 - cos(dPhi) );
    
}

double BranchingBendingCosine::energy(Bead* b1, Bead* b2, Bead* b3, Bead* b4, double k_bend, double phi0, double d){
    
    double L1 = sqrt(ScalarProductStretched(b1->coordinate, b1->force, b2->coordinate, b2->force,
                                            b1->coordinate, b1->force, b2->coordinate, b2->force, d));
    double L2 = sqrt(ScalarProductStretched(b3->coordinate, b3->force, b4->coordinate, b4->force,
                                            b3->coordinate, b3->force, b4->coordinate, b4->force, d));
    
    double L1L2 = L1*L2;
    double l1l2 = ScalarProductStretched(b1->coordinate, b1->force, b2->coordinate, b2->force,
                                         b3->coordinate, b3->force, b4->coordinate, b4->force, d);
    
    double phi = acos(l1l2 / L1L2);
    double dPhi = phi-phi0;
    
    return k_bend * ( 1 - cos(dPhi) );
    
}

void BranchingBendingCosine::forces(Bead* b1, Bead* b2, Bead* b3, Bead* b4, double k_bend, double phi0 ){
    
    //l1 = b2-b3; l2 = b4-b3;
    
    double L1 = sqrt(ScalarProduct(b1->coordinate, b2->coordinate, b1->coordinate, b2->coordinate));
    double L2 = sqrt(ScalarProduct(b3->coordinate, b4->coordinate, b3->coordinate, b4->coordinate));
    double l1l2 = ScalarProduct(b1->coordinate, b2->coordinate, b3->coordinate, b4->coordinate);
    
    //invL = 1/L;
    double invL1 = 1/L1;
    double invL2 = 1/L2;
    double A = invL1*invL2;
    double B = l1l2*invL1*A*A*L2;
    double C = l1l2*invL2*A*A*L1;
    
    double phi = acos(l1l2 *A);
    double dPhi = phi-phi0;
    
    double k =  k_bend* sin(dPhi)/sin(phi);
    
    //force on i, f = k*(-A*l2 + 2*B*l1):
    b1->force[0] +=  k * ( (b3->coordinate[0] - b4->coordinate[0])*A + (b2->coordinate[0] - b1->coordinate[0])*B );
    b1->force[1] +=  k * ( (b3->coordinate[1] - b4->coordinate[1])*A + (b2->coordinate[1] - b1->coordinate[1])*B );
    b1->force[2] +=  k * ( (b3->coordinate[2] - b4->coordinate[2])*A + (b2->coordinate[2] - b1->coordinate[2])*B );
    
    
    //force on i+1, f = k*(A*l2 - 2*B*l1):
    b2->force[0] +=  k * ( (-b3->coordinate[0] + b4->coordinate[0])*A - (b2->coordinate[0] - b1->coordinate[0])*B );
    b2->force[1] +=  k * ( (-b3->coordinate[1] + b4->coordinate[1])*A - (b2->coordinate[1] - b1->coordinate[1])*B );
    b2->force[2] +=  k * ( (-b3->coordinate[2] + b4->coordinate[2])*A - (b2->coordinate[2] - b1->coordinate[2])*B );
    
    //force on j, k*(-A*l1 + 2*C*l2):
    b3->force[0] +=  k *( (b1->coordinate[0] - b2->coordinate[0])*A + (b4->coordinate[0] - b3->coordinate[0])*C );
    b3->force[1] +=  k *( (b1->coordinate[1] - b2->coordinate[1])*A + (b4->coordinate[1] - b3->coordinate[1])*C );
    b3->force[2] +=  k *( (b1->coordinate[2] - b2->coordinate[2])*A + (b4->coordinate[2] - b3->coordinate[2])*C );
    
    //force on j+1, k*(A*l1 - 2*C*l2):
    b4->force[0] +=  k *( (-b1->coordinate[0] + b2->coordinate[0])*A - (b4->coordinate[0] - b3->coordinate[0])*C );
    b4->force[1] +=  k *( (-b1->coordinate[1] + b2->coordinate[1])*A - (b4->coordinate[1] - b3->coordinate[1])*C );
    b4->force[2] +=  k *( (-b1->coordinate[2] + b2->coordinate[2])*A - (b4->coordinate[2] - b3->coordinate[2])*C );
    
    
}

void BranchingBendingCosine::forcesAux(Bead* b1, Bead* b2, Bead* b3, Bead* b4, double k_bend, double phi0 ){
    
    //l1 = b2-b3; l2 = b4-b3;
    
    double L1 = sqrt(ScalarProduct(b1->coordinateAux, b2->coordinateAux, b1->coordinateAux, b2->coordinateAux));
    double L2 = sqrt(ScalarProduct(b3->coordinateAux, b4->coordinateAux, b3->coordinateAux, b4->coordinateAux));
    double l1l2 = ScalarProduct(b1->coordinateAux, b2->coordinateAux, b3->coordinateAux, b4->coordinateAux);
    
    //invL = 1/L;
    double invL1 = 1/L1;
    double invL2 = 1/L2;
    double A = invL1*invL2;
    double B = l1l2*invL1*A*A*L2;
    double C = l1l2*invL2*A*A*L1;
    
    double phi = acos(l1l2 *A);
    double dPhi = phi-phi0;
    
    double k =  k_bend* sin(dPhi)/sin(phi);
    
    //force on i, f = k*(-A*l2 + 2*B*l1):
    b1->forceAux[0] +=  k * ( (b3->coordinateAux[0] - b4->coordinateAux[0])*A + (b2->coordinateAux[0] - b1->coordinateAux[0])*B );
    b1->forceAux[1] +=  k * ( (b3->coordinate[1] - b4->coordinate[1])*A + (b2->coordinate[1] - b1->coordinateAux[1])*B );
    b1->forceAux[2] +=  k * ( (b3->coordinateAux[2] - b4->coordinateAux[2])*A + (b2->coordinateAux[2] - b1->coordinateAux[2])*B );
    
    
    //force on i+1, f = k*(A*l2 - 2*B*l1):
    b2->forceAux[0] +=  k * ( (-b3->coordinateAux[0] + b4->coordinateAux[0])*A - (b2->coordinateAux[0] - b1->coordinateAux[0])*B );
    b2->forceAux[1] +=  k * ( (-b3->coordinateAux[1] + b4->coordinateAux[1])*A - (b2->coordinateAux[1] - b1->coordinateAux[1])*B );
    b2->forceAux[2] +=  k * ( (-b3->coordinateAux[2] + b4->coordinateAux[2])*A - (b2->coordinateAux[2] - b1->coordinateAux[2])*B );
    
    //force on j, k*(-A*l1 + 2*C*l2):
    b3->forceAux[0] +=  k *( (b1->coordinateAux[0] - b2->coordinateAux[0])*A + (b4->coordinateAux[0] - b3->coordinateAux[0])*C );
    b3->forceAux[1] +=  k *( (b1->coordinateAux[1] - b2->coordinateAux[1])*A + (b4->coordinateAux[1] - b3->coordinateAux[1])*C );
    b3->forceAux[2] +=  k *( (b1->coordinateAux[2] - b2->coordinateAux[2])*A + (b4->coordinateAux[2] - b3->coordinateAux[2])*C );
    
    //force on j+1, k*(A*l1 - 2*C*l2):
    b4->forceAux[0] +=  k *( (-b1->coordinateAux[0] + b2->coordinateAux[0])*A - (b4->coordinateAux[0] - b3->coordinateAux[0])*C );
    b4->forceAux[1] +=  k *( (-b1->coordinateAux[1] + b2->coordinateAux[1])*A - (b4->coordinateAux[1] - b3->coordinateAux[1])*C );
    b4->forceAux[2] +=  k *( (-b1->coordinateAux[2] + b2->coordinateAux[2])*A - (b4->coordinateAux[2] - b3->coordinateAux[2])*C );
    
}