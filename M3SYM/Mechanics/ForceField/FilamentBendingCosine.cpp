
//------------------------------------------------------------------
//  **M3SYM** - Simulation Package for the Mechanochemical
//              Dynamics of Active Networks, 3rd Generation
//
//  Copyright (2014) Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the Papoian lab page for installation and documentation:
//  http://papoian.chem.umd.edu/
//------------------------------------------------------------------

#include "FilamentBendingCosine.h"

#include "MathFunctions.h"
#include "Bead.h"

using namespace mathfunc;

double FilamentBendingCosine::energy(Bead* b1, Bead* b2, Bead* b3,
                                     double kBend, double eqTheta){
    
    double L1 = sqrt(scalarProduct(b1->coordinate, b2->coordinate,
                                   b1->coordinate, b2->coordinate));
    double L2 = sqrt(scalarProduct(b2->coordinate, b3->coordinate,
                                   b2->coordinate, b3->coordinate));
    
    double L1L2 = L1*L2;
    double l1l2 = scalarProduct(b1->coordinate, b2->coordinate,
                                b2->coordinate, b3->coordinate);
    
    double phi = acos(l1l2 / L1L2);
    double dPhi = phi-eqTheta;
    
    return kBend * ( 1 - cos(dPhi) );
}

double FilamentBendingCosine::energy(Bead* b1, Bead* b2, Bead* b3,
                                     double kBend, double eqTheta, double d ){
    
    double L1 = sqrt(scalarProductStretched(b1->coordinate, b1->force,
                                            b2->coordinate, b2->force,
                                            b1->coordinate, b1->force,
                                            b2->coordinate, b2->force, d));
    double L2 = sqrt(scalarProductStretched(b2->coordinate, b2->force,
                                            b3->coordinate, b3->force,
                                            b2->coordinate, b2->force,
                                            b3->coordinate, b3->force, d));
    
    double L1L2 = L1*L2;
    double l1l2 = scalarProductStretched(b1->coordinate, b1->force,
                                         b2->coordinate, b2->force,
                                         b2->coordinate, b2->force,
                                         b3->coordinate, b3->force, d);
    
    double phi = acos(l1l2 / L1L2);
    double dPhi = phi-eqTheta;
    
    return kBend * ( 1 - cos(dPhi) );
}

void FilamentBendingCosine::forces(Bead* b1, Bead* b2, Bead* b3,
                                   double kBend, double eqTheta ){
    
    double k = 0;
    double L1 = sqrt(scalarProduct(b1->coordinate, b2->coordinate,
                                   b1->coordinate, b2->coordinate));
    double L2 = sqrt(scalarProduct(b2->coordinate, b3->coordinate,
                                   b2->coordinate, b3->coordinate));
    double l1l2 = scalarProduct(b1->coordinate, b2->coordinate,
                                b2->coordinate, b3->coordinate);
    
    //invL = 1/L;
    double invL1 = 1/L1;
    double invL2 = 1/L2;
    double A = invL1*invL2;
    double B = l1l2*invL1*A*A*L2;
    double C = l1l2*invL2*A*A*L1;
   
    if (eqTheta == 0) {k = kBend;}
    
    else{
        double phi = acos(l1l2 *A);
        double dPhi = phi-eqTheta;
        
        k =  kBend* sin(dPhi)/sin(phi);
    }
    
    //force on i-1, f = k*(-A*l2 + B*l1):
    b1->force[0] +=  k * ((-b3->coordinate[0] + b2->coordinate[0])*A +
                          (b2->coordinate[0] - b1->coordinate[0])*B );
    b1->force[1] +=  k * ((-b3->coordinate[1] + b2->coordinate[1])*A +
                          (b2->coordinate[1] - b1->coordinate[1])*B );
    b1->force[2] +=  k * ((-b3->coordinate[2] + b2->coordinate[2])*A +
                          (b2->coordinate[2] - b1->coordinate[2])*B );
    
    
    //force on i, f = k*(A*(l1-l2) - B*l1 + C*l2):
    b2->force[0] +=  k *( (b3->coordinate[0] -
                           2*b2->coordinate[0] +
                           b1->coordinate[0])*A -
                           (b2->coordinate[0] -
                           b1->coordinate[0])*B +
                           (b3->coordinate[0] -
                           b2->coordinate[0])*C );
    
    b2->force[1] +=  k *( (b3->coordinate[1] -
                           2*b2->coordinate[1] +
                           b1->coordinate[1])*A -
                           (b2->coordinate[1] -
                           b1->coordinate[1])*B +
                           (b3->coordinate[1] -
                           b2->coordinate[1])*C );
    
    b2->force[2] +=  k *( (b3->coordinate[2] -
                           2*b2->coordinate[2] +
                           b1->coordinate[2])*A -
                           (b2->coordinate[2] -
                           b1->coordinate[2])*B +
                           (b3->coordinate[2] -
                           b2->coordinate[2])*C );
    
    //force on i-1, f = k*(A*l - B*l2):
    b3->force[0] +=  k *( (b2->coordinate[0] -
                           b1->coordinate[0])*A -
                           (b3->coordinate[0] -
                           b2->coordinate[0])*C );
    
    b3->force[1] +=  k *( (b2->coordinate[1] -
                           b1->coordinate[1])*A -
                           (b3->coordinate[1] -
                           b2->coordinate[1])*C );
    
    b3->force[2] +=  k *( (b2->coordinate[2] -
                           b1->coordinate[2])*A -
                           (b3->coordinate[2] -
                           b2->coordinate[2])*C );
}

void FilamentBendingCosine::forcesAux(Bead* b1, Bead* b2, Bead* b3,
                                      double kBend, double eqTheta){
    
    double k = 0;
    double L1 = sqrt(scalarProduct(b1->coordinateAux, b2->coordinateAux,
                                   b1->coordinateAux, b2->coordinateAux));
    double L2 = sqrt(scalarProduct(b2->coordinateAux, b3->coordinateAux,
                                   b2->coordinateAux, b3->coordinateAux));
    double l1l2 = scalarProduct(b1->coordinateAux, b2->coordinateAux,
                                b2->coordinateAux, b3->coordinateAux);
    
    //invL = 1/L;
    double invL1 = 1/L1;
    double invL2 = 1/L2;
    double A = invL1*invL2;
    double B = l1l2*invL1*A*A*L2;
    double C = l1l2*invL2*A*A*L1;
    
    if (eqTheta == 0) {k = kBend;}
    
    else{
        double phi = acos(l1l2 *A);
        double dPhi = phi-eqTheta;
        
        k =  kBend* sin(dPhi)/sin(phi);
    }
    
    //force on i-1, f = k*(-A*l2 + B*l1):
    b1->force[0] +=  k * ( (-b3->coordinateAux[0] +
                            b2->coordinateAux[0])*A +
                            (b2->coordinateAux[0] -
                            b1->coordinateAux[0])*B );
    
    b1->force[1] +=  k * ( (-b3->coordinateAux[1] +
                            b2->coordinateAux[1])*A +
                            (b2->coordinateAux[1] -
                            b1->coordinateAux[1])*B );
    
    b1->force[2] +=  k* ( (-b3->coordinateAux[2] +
                           b2->coordinateAux[2])*A +
                           (b2->coordinateAux[2] -
                           b1->coordinateAux[2])*B );
    
    //force on i, f = k*(A*(l1-l2) - B*l1 + C*l2):
    b2->forceAux[0] +=  k *( (b3->coordinateAux[0] -
                              2*b2->coordinateAux[0] +
                              b1->coordinateAux[0])*A -
                              (b2->coordinateAux[0] -
                              b1->coordinateAux[0])*B +
                              (b3->coordinateAux[0] -
                              b2->coordinateAux[0])*C );
    
    b2->forceAux[1] +=  k *( (b3->coordinateAux[1] -
                              2*b2->coordinateAux[1] +
                              b1->coordinateAux[1])*A -
                              (b2->coordinateAux[1] -
                              b1->coordinateAux[1])*B +
                              (b3->coordinateAux[1] -
                              b2->coordinateAux[1])*C );
    
    b2->forceAux[2] +=  k *( (b3->coordinateAux[2] -
                              2*b2->coordinateAux[2] +
                              b1->coordinateAux[2])*A -
                              (b2->coordinateAux[2] -
                              b1->coordinateAux[2])*B +
                              (b3->coordinateAux[2] -
                              b2->coordinateAux[2])*C );
    
    //force on i-1, f = k*(A*l - B*l2):
    b3->forceAux[0] +=  k *( (b2->coordinateAux[0] -
                              b1->coordinateAux[0])*A -
                              (b3->coordinateAux[0] -
                              b2->coordinateAux[0])*C );
    
    b3->forceAux[1] +=  k *( (b2->coordinateAux[1] -
                              b1->coordinateAux[1])*A -
                              (b3->coordinateAux[1] -
                              b2->coordinateAux[1])*C );
    
    b3->forceAux[2] +=  k *( (b2->coordinateAux[2] -
                              b1->coordinateAux[2])*A -
                              (b3->coordinateAux[2] -
                              b2->coordinateAux[2])*C );
    
}