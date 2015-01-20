
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

#include <cmath>

#include "BranchingPositionCosine.h"

#include "Bead.h"

#include "MathFunctions.h"

////!!!!!!!!!!NOT READY THIS IS JUST A TEMPLATE!!!!!!! BUT SHOULD COMPILE
using namespace mathfunc;

double BranchingPositionCosine::energy(Bead* b1, Bead* b2, Bead* b3, Bead* b4,
                                       double kBend, double eqTheta){
    
    double L1 = sqrt(scalarProduct(b1->coordinate, b2->coordinate,
                                   b1->coordinate, b2->coordinate));
    double L2 = sqrt(scalarProduct(b3->coordinate, b4->coordinate,
                                   b3->coordinate, b4->coordinate));
    
    double L1L2 = L1*L2;
    double l1l2 = scalarProduct(b1->coordinate, b2->coordinate,
                                b3->coordinate, b4->coordinate);
    
    double theta = acos(l1l2 / L1L2);
    double dtheta = theta-eqTheta;
    
    double U = kBend * ( 1 - cos(dtheta) );
    
    return 0;
    
}

double BranchingPositionCosine::energy(Bead* b1, Bead* b2, Bead* b3, Bead* b4,
                                       double kBend, double eqTheta, double d){
    
    double L1 = sqrt(scalarProductStretched(b1->coordinate, b1->force,
                                            b2->coordinate, b2->force,
                                            b1->coordinate, b1->force,
                                            b2->coordinate, b2->force, d));
    double L2 = sqrt(scalarProductStretched(b3->coordinate, b3->force,
                                            b4->coordinate, b4->force,
                                            b3->coordinate, b3->force,
                                            b4->coordinate, b4->force, d));
    
    double L1L2 = L1*L2;
    double l1l2 = scalarProductStretched(b1->coordinate, b1->force,
                                         b2->coordinate, b2->force,
                                         b3->coordinate, b3->force,
                                         b4->coordinate, b4->force, d);
    
    double theta = acos(l1l2 / L1L2);
    double dtheta = theta-eqTheta;
    double U = kBend * ( 1 - cos(dtheta) );
    
    return 0;
    
}

void BranchingPositionCosine::forces(Bead* b1, Bead* b2, Bead* b3, Bead* b4,
                                     double kBend, double eqTheta ){
    
    //l1 = b2-b3; l2 = b4-b3;
    
    double L1 = sqrt(scalarProduct(b1->coordinate, b2->coordinate,
                                   b1->coordinate, b2->coordinate));
    double L2 = sqrt(scalarProduct(b3->coordinate, b4->coordinate,
                                   b3->coordinate, b4->coordinate));
    double l1l2 = scalarProduct(b1->coordinate, b2->coordinate,
                                b3->coordinate, b4->coordinate);
    
    //invL = 1/L;
    double invL1 = 1/L1;
    double invL2 = 1/L2;
    double A = invL1*invL2;
    double B = l1l2*invL1*A*A*L2;
    double C = l1l2*invL2*A*A*L1;
    
    double theta = acos(l1l2 *A);
    double dtheta = theta-eqTheta;
    
    double k =  kBend* sin(dtheta)/sin(theta);
    
    //force on i, f = k*(-A*l2 + 2*B*l1):
    
    
    //force on i+1, f = k*(A*l2 - 2*B*l1):

    
    //force on j, k*(-A*l1 + 2*C*l2):
    
    //force on j+1, k*(A*l1 - 2*C*l2):
    
}

void BranchingPositionCosine::forcesAux(Bead* b1, Bead* b2, Bead* b3, Bead* b4,
                                        double kBend, double eqTheta ){
    
    //l1 = b2-b3; l2 = b4-b3;
    
    double L1 = sqrt(scalarProduct(b1->coordinateAux, b2->coordinateAux,
                                   b1->coordinateAux, b2->coordinateAux));
    double L2 = sqrt(scalarProduct(b3->coordinateAux, b4->coordinateAux,
                                   b3->coordinateAux, b4->coordinateAux));
    double l1l2 = scalarProduct(b1->coordinateAux, b2->coordinateAux,
                                b3->coordinateAux, b4->coordinateAux);
    
    //invL = 1/L;
    double invL1 = 1/L1;
    double invL2 = 1/L2;
    double A = invL1*invL2;
    double B = l1l2*invL1*A*A*L2;
    double C = l1l2*invL2*A*A*L1;
    
    double theta = acos(l1l2 *A);
    double dtheta = theta-eqTheta;
    
    double k =  kBend* sin(dtheta)/sin(theta);
    
    //force on i, f = k*(-A*l2 + 2*B*l1):
    
    
    //force on i+1, f = k*(A*l2 - 2*B*l1):
   
    
    //force on j, k*(-A*l1 + 2*C*l2):
   
    
    //force on j+1, k*(A*l1 - 2*C*l2):
       
}