
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v3.1
//
//  Copyright (2015-2016)  Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the MEDYAN web page for more information:
//  http://www.medyan.org
//------------------------------------------------------------------

#include "FilamentBendingHarmonic.h"

#include "Bead.h"

#include "MathFunctions.h"

using namespace mathfunc;
#ifdef CROSSCHECK
double FilamentBendingHarmonic::energy(Bead* b1, Bead* b2, Bead* b3,
                                       double kBend, double eqTheta){
    
    double L1 = sqrt(scalarProduct(b1->coordinate, b2->coordinate,
                                   b1->coordinate, b2->coordinate));
    double L2 = sqrt(scalarProduct(b2->coordinate, b3->coordinate,
                                   b2->coordinate, b3->coordinate));
    
    double L1L2 = L1*L2;
    double l1l2 = scalarProduct(b1->coordinate, b2->coordinate,
                                b2->coordinate, b3->coordinate);

    return kBend * ( 1 - l1l2 / L1L2 );
    
}

double FilamentBendingHarmonic::energy(Bead* b1, Bead* b2, Bead* b3,
                                       double kBend, double eqTheta, double d){
    
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

    return kBend * ( 1 - l1l2 / L1L2 );
    
}

void FilamentBendingHarmonic::forces(Bead* b1, Bead* b2, Bead* b3,
                                     double kBend, double eqTheta){
    
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
    
    //force on i-1, f = k*(-A*l2 + B*l1):
    b1->force[0] +=  kBend * ( (-b3->coordinate[0] + b2->coordinate[0])*A +
                              (b2->coordinate[0] - b1->coordinate[0])*B );
    b1->force[1] +=  kBend * ( (-b3->coordinate[1] + b2->coordinate[1])*A +
                              (b2->coordinate[1] - b1->coordinate[1])*B );
    b1->force[2] +=  kBend * ( (-b3->coordinate[2] + b2->coordinate[2])*A +
                              (b2->coordinate[2] - b1->coordinate[2])*B );
    
    
    //force on i, f = k*(A*(l1-l2) - B*l1 + C*l2):
    b2->force[0] +=  kBend *( (b3->coordinate[0] - 2*b2->coordinate[0] +
                               b1->coordinate[0])*A -
                             (b2->coordinate[0] - b1->coordinate[0])*B +
                             (b3->coordinate[0] - b2->coordinate[0])*C );
    
    b2->force[1] +=  kBend *( (b3->coordinate[1] - 2*b2->coordinate[1] +
                               b1->coordinate[1])*A -
                             (b2->coordinate[1] - b1->coordinate[1])*B +
                             (b3->coordinate[1] - b2->coordinate[1])*C );
    
    b2->force[2] +=  kBend *( (b3->coordinate[2] - 2*b2->coordinate[2] +
                               b1->coordinate[2])*A -
                             (b2->coordinate[2] - b1->coordinate[2])*B +
                             (b3->coordinate[2] - b2->coordinate[2])*C );
    
    //force on i-1, f = k*(A*l - B*l2):
    b3->force[0] +=  kBend *( (b2->coordinate[0] - b1->coordinate[0])*A -
                             (b3->coordinate[0] - b2->coordinate[0])*C );
    b3->force[1] +=  kBend *( (b2->coordinate[1] - b1->coordinate[1])*A -
                             (b3->coordinate[1] - b2->coordinate[1])*C );
    b3->force[2] +=  kBend *( (b2->coordinate[2] - b1->coordinate[2])*A -
                             (b3->coordinate[2] - b2->coordinate[2])*C );
    
    
}

void FilamentBendingHarmonic::forcesAux(Bead* b1, Bead* b2, Bead* b3,
                                        double kBend, double eqTheta){
    
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
    
    //force on i-1, f = k*(-A*l2 + B*l1):
    b1->forceAux[0] +=  kBend * ( (-b3->coordinate[0] + b2->coordinate[0])*A +
                                 (b2->coordinate[0] - b1->coordinate[0])*B );
    b1->forceAux[1] +=  kBend * ( (-b3->coordinate[1] + b2->coordinate[1])*A +
                                 ( b2->coordinate[1] - b1->coordinate[1])*B );
    b1->forceAux[2] +=  kBend * ( (-b3->coordinate[2] + b2->coordinate[2])*A +
                                 (b2->coordinate[2] - b1->coordinate[2])*B );
    
    //force on i, f = k*(A*(l1-l2) - B*l1 + C*l2):
    b2->forceAux[0] +=  kBend *( (b3->coordinate[0] - 2*b2->coordinate[0] +
                                  b1->coordinate[0])*A -
                                (b2->coordinate[0] - b1->coordinate[0])*B +
                                (b3->coordinate[0] - b2->coordinate[0])*C );
    
    b2->forceAux[1] +=  kBend *( (b3->coordinate[1] - 2*b2->coordinate[1] +
                                  b1->coordinate[1])*A -
                                (b2->coordinate[1] - b1->coordinate[1])*B +
                                (b3->coordinate[1] - b2->coordinate[1])*C );
    
    b2->forceAux[2] +=  kBend *( (b3->coordinate[2] - 2*b2->coordinate[2] +
                                  b1->coordinate[2])*A -
                                (b2->coordinate[2] - b1->coordinate[2])*B +
                                (b3->coordinate[2] - b2->coordinate[2])*C );
    
    //force on i-1, f = k*(A*l - B*l2):
    b3->forceAux[0] +=  kBend *( (b2->coordinate[0] - b1->coordinate[0])*A -
                                (b3->coordinate[0] - b2->coordinate[0])*C );
    b3->forceAux[1] +=  kBend *( (b2->coordinate[1] - b1->coordinate[1])*A -
                                (b3->coordinate[1] - b2->coordinate[1])*C );
    b3->forceAux[2] +=  kBend *( (b2->coordinate[2] - b1->coordinate[2])*A -
                                (b3->coordinate[2] - b2->coordinate[2])*C );
    
}
#endif;
