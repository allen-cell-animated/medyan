
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

#include <cmath>

#include "BranchingBendingCosine.h"
#include "BranchingBending.h"

#include "BranchingPoint.h"
#include "Bead.h"

#include "MathFunctions.h"

using namespace mathfunc;

double BranchingBendingCosine::energy(double *coord, double *f, int *beadSet,
                                      double *kbend, double *eqt){

    int n = BranchingBending<BranchingBendingCosine>::n;
    int nint = BranchingPoint::getBranchingPoints().size();
    
    double *coord1, *coord2, *coord3, *coord4, dist, U_i, L1, L2, L1L2, l1l2, phi, dPhi;
    
    double U = 0;
    
    for(int i = 0; i < nint; i += 1) {
        
        coord1 = &coord[3 * beadSet[n * i]];
        coord2 = &coord[3 * beadSet[n * i + 1]];
        coord3 = &coord[3 * beadSet[n * i + 2]];
        coord4 = &coord[3 * beadSet[n * i + 3]];
        
        L1 = sqrt(scalarProduct(coord1, coord2,
                                coord1, coord2));
        L2 = sqrt(scalarProduct(coord3, coord4,
                                coord3, coord4));
        
        L1L2 = L1*L2;
        l1l2 = scalarProduct(coord1, coord2,
                             coord3, coord4);
        
        phi = safeacos(l1l2 / L1L2);
        dPhi = phi-eqt[i];
        
        U_i = kbend[i] * ( 1 - cos(dPhi) );
        
        if(fabs(U_i) == numeric_limits<double>::infinity()
           || U_i != U_i || U_i < -1.0) {
            
            //set culprit and return
            BranchingInteractions::_branchingCulprit = BranchingPoint::getBranchingPoints()[i];
            
            return -1;
        }
        
        U += U_i;
    }
    
    return U;
}

double BranchingBendingCosine::energy(double *coord, double *f, int *beadSet,
                                      double *kbend, double *eqt, double d){
    
    int n = BranchingBending<BranchingBendingCosine>::n;
    int nint = BranchingPoint::getBranchingPoints().size();
    
    double *coord1, *coord2, *coord3, *coord4, *force1, *force2, *force3, *force4, dist, U_i, L1, L2, L1L2, l1l2, phi, dPhi;
    
    double U = 0;
    
    for(int i = 0; i < nint; i += 1) {
        
        coord1 = &coord[3 * beadSet[n * i]];
        coord2 = &coord[3 * beadSet[n * i + 1]];
        coord3 = &coord[3 * beadSet[n * i + 2]];
        coord3 = &coord[3 * beadSet[n * i + 3]];
        
        force1 = &f[3 * beadSet[n * i]];
        force2 = &f[3 * beadSet[n * i + 1]];
        force3 = &f[3 * beadSet[n * i + 2]];
        force4 = &f[3 * beadSet[n * i + 3]];
        
        
        L1 = sqrt(scalarProductStretched(coord1, force1, coord2, force2,
                                         coord1, force1, coord2, force2, d));
        L2 = sqrt(scalarProductStretched(coord3, force3, coord4, force4,
                                         coord3, force3, coord4, force4, d));
        
        L1L2 = L1*L2;
        l1l2 = scalarProductStretched(coord1, force1, coord2, force2,
                                      coord3, force3, coord4, force4, d);
        
        phi = safeacos(l1l2 / L1L2);
        dPhi = phi-eqt[i];
        
        U_i = kbend[i] * ( 1 - cos(dPhi) );
        
        if(fabs(U_i) == numeric_limits<double>::infinity()
           || U_i != U_i || U_i < -1.0) {
            
            //set culprit and return
            BranchingInteractions::_branchingCulprit = BranchingPoint::getBranchingPoints()[i];
            
            return -1;
        }
        
        U += U_i;
    }
    
    return U;
}

void BranchingBendingCosine::forces(double *coord, double *f, int *beadSet,
                                    double *kbend, double *eqt){
    
   
    int n = BranchingBending<BranchingBendingCosine>::n;
    int nint = BranchingPoint::getBranchingPoints().size();
    
    double *coord1, *coord2, *coord3, *coord4, *force1, *force2, *force3, *force4;
    double dist, U_i, L1, L2, L1L2, l1l2, phi, dPhi, A, B, C, invL1, invL2, k;
    
    double U = 0;
    
    for(int i = 0; i < nint; i += 1) {
        
        coord1 = &coord[3 * beadSet[n * i]];
        coord2 = &coord[3 * beadSet[n * i + 1]];
        coord3 = &coord[3 * beadSet[n * i + 2]];
        coord3 = &coord[3 * beadSet[n * i + 3]];
        
        force1 = &f[3 * beadSet[n * i]];
        force2 = &f[3 * beadSet[n * i + 1]];
        force3 = &f[3 * beadSet[n * i + 2]];
        force4 = &f[3 * beadSet[n * i + 3]];
        
        
        L1 = sqrt(scalarProduct(coord1, coord2,
                                coord1, coord2));
        L2 = sqrt(scalarProduct(coord3, coord4,
                                coord3, coord4));
        
        L1L2 = L1*L2;
        l1l2 = scalarProduct(coord1, coord2,
                             coord3, coord4);
        
        invL1 = 1/L1;
        invL2 = 1/L2;
        A = invL1*invL2;
        B = l1l2*invL1*A*A*L2;
        C = l1l2*invL2*A*A*L1;
        
        phi = safeacos(l1l2 / L1L2);
        dPhi = phi-eqt[i];
        
        k =  kbend[i] * sin(dPhi)/sin(phi);
        
        //force on i, f = k*(-A*l2 + 2*B*l1):
        force1[0] += k * ((coord3[0] - coord4[0])*A +
                          (coord2[0] - coord1[0])*B );
        force1[1] += k * ((coord3[1] - coord4[1])*A +
                          (coord2[1] - coord1[1])*B );
        force1[2] += k * ((coord3[2] - coord4[2])*A +
                          (coord2[2] - coord1[2])*B );
        
        
        //force on i+1, f = k*(A*l2 - 2*B*l1):
        force2[0] += k * ((-coord3[0] + coord4[0])*A -
                          (coord2[0] - coord1[0])*B );
        force2[1] += k * ((-coord3[1] + coord4[1])*A -
                          (coord2[1] - coord1[1])*B );
        force2[2] += k * ((-coord3[2] + coord4[2])*A -
                          (coord2[2] - coord1[2])*B );
        
        //force on j, k*(-A*l1 + 2*C*l2):
        force3[0] += k *((coord1[0] - coord2[0])*A +
                         (coord4[0] - coord3[0])*C );
        force3[1] += k *((coord1[1] - coord2[1])*A +
                         (coord4[1] - coord3[1])*C );
        force3[2] += k *((coord1[2] - coord2[2])*A +
                         (coord4[2] - coord3[2])*C );
        
        //force on j+1, k*(A*l1 - 2*C*l2):
        force4[0] += k *((-coord1[0] + coord2[0])*A -
                         (coord4[0] - coord3[0])*C );
        force4[1] += k *((-coord1[1] + coord2[1])*A -
                         (coord4[1] - coord3[1])*C );
        force4[2] += k *((-coord1[2] + coord2[2])*A -
                         (coord4[2] - coord3[2])*C );
    }
}
