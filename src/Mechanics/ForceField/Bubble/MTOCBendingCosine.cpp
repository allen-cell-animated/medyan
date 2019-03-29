
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

#include "MTOCBendingCosine.h"
#include "MTOCBending.h"
#include "Bead.h"
#include "Bubble.h"
#include "MTOC.h"

#include "MathFunctions.h"

using namespace mathfunc;

double MTOCBendingCosine::energy(double *coord, double *f, int *beadSet,
                                      double *kbend, double radius){
    
    double *coord1, *coord2, *coord3, U_i, L1, L2, L1L2, l1l2, phi, dPhi;
    double U = 0.0;
    
    int n = MTOCBending<MTOCBendingCosine>::n;
    for(auto mtoc : MTOC::getMTOCs()) {
        int nint = mtoc->getFilaments().size();
        
        coord1 = &coord[3 * beadSet[0]]; //coordinate of MTOC
        
        for(int i = 0; i < nint; i+=1){
            coord2 = &coord[3 * beadSet[n * i + 1]];
            coord3 = &coord[3 * beadSet[n * i + 2]];
            
            L1 = sqrt(scalarProduct(coord1, coord2,
                                    coord1, coord2));
            L2 = sqrt(scalarProduct(coord2, coord3,
                                    coord2, coord3));
            
            L1L2 = L1*L2;
            l1l2 = scalarProduct(coord1, coord2,
                                 coord2, coord3);
            //the equilibrium theta is 0
            dPhi = safeacos(l1l2 / L1L2);
            
            U_i = kbend[i] * ( 1 - cos(dPhi) );
            
            if(fabs(U_i) == numeric_limits<double>::infinity()
               || U_i != U_i || U_i < -1.0) {
                
                //set culprit and return TODO
                //FilamentInteractions::_filamentCulprit = Filament::getFilaments()[i];
                
                return -1;
            }
            
            U += U_i;
            
        }
    }
    return U;
    
}

double MTOCBendingCosine::energy(double *coord, double *f, int *beadSet,
                                      double *kbend, double radius, double d){
    
    //Do not use for now
    
    double *coord1, *coord2, *coord3, U_i, L1, L2, L1L2, l1l2, phi, dPhi;
    double U = 0.0;
    
    int n = MTOCBending<MTOCBendingCosine>::n;
    for(auto mtoc : MTOC::getMTOCs()) {
        int nint = mtoc->getFilaments().size();
        
        coord1 = &coord[3 * beadSet[0]]; //coordinate of MTOC
        
        for(int i = 0; i < nint; i+=1){
            coord2 = &coord[3 * beadSet[n * i + 1]];
            coord3 = &coord[3 * beadSet[n * i + 2]];
            
            
            L1 = sqrt(scalarProduct(coord1, coord2,
                                    coord1, coord2));
            L2 = sqrt(scalarProduct(coord2, coord3,
                                    coord2, coord3));
            
            L1L2 = L1*L2;
            l1l2 = scalarProduct(coord1, coord2,
                                 coord2, coord3);
            //the equilibrium theta is 0
            dPhi = safeacos(l1l2 / L1L2);
            
            U_i = kbend[i] * ( 1 - cos(dPhi) );
            
            if(fabs(U_i) == numeric_limits<double>::infinity()
               || U_i != U_i || U_i < -1.0) {
                
                //set culprit and return TODO
                //FilamentInteractions::_filamentCulprit = Filament::getFilaments()[i];
                
                return -1;
            }
            
            U += U_i;
            
        }
    }
    return U;
    
}

void MTOCBendingCosine::forces(double *coord, double *f, int *beadSet,
                                    double *kbend, double radius){
    
    
    int n = MTOCBending<MTOCBendingCosine>::n;
    for(auto mtoc : MTOC::getMTOCs()) {
        int nint = mtoc->getFilaments().size();
        
        double *coord1, *coord2, *coord3, *force1, *force2, *force3,
        L1, L2, l1l2, invL1, invL2, A,B,C, phi, dPhi, k;
        
        coord1 = &coord[3 * beadSet[0]]; //coordinate of MTOC
        force1 = &f[3 * beadSet[0]];
        
        for(int i = 0; i < nint; i+=1){
            coord2 = &coord[3 * beadSet[n * i + 1]];
            coord3 = &coord[3 * beadSet[n * i + 2]];
            
            force2 = &f[3 * beadSet[n * i + 1]];
            force3 = &f[3 * beadSet[n * i + 2]];
            
            L1 = sqrt(scalarProduct(coord1, coord2,
                                    coord1, coord2));
            L2 = sqrt(scalarProduct(coord2, coord3,
                                    coord2, coord3));
            
            l1l2 = scalarProduct(coord1, coord2,
                                 coord2, coord3);
            
            invL1 = 1/L1;
            invL2 = 1/L2;
            A = invL1*invL2;
            B = l1l2*invL1*A*A*L2;
            C = l1l2*invL2*A*A*L1;
            
            k = kbend[i];
            
            force1[0] +=  k * ((-coord3[0] + coord2[0])*A +
                               (coord2[0] - coord1[0])*B );
            force1[1] +=  k * ((-coord3[1] + coord2[1])*A +
                               (coord2[1] - coord1[1])*B );
            force1[2] +=  k * ((-coord3[2] + coord2[2])*A +
                               (coord2[2] - coord1[2])*B );
            
            
            //force on i, f = k*(A*(l1-l2) - B*l1 + C*l2):
            force2[0] +=  k *( (coord3[0] - 2*coord2[0] + coord1[0])*A -
                              (coord2[0] - coord1[0])*B +
                              (coord3[0] - coord2[0])*C );
            
            force2[1] +=  k *( (coord3[1] - 2*coord2[1] + coord1[1])*A -
                              (coord2[1] - coord1[1])*B +
                              (coord3[1] - coord2[1])*C );
            
            force2[2] +=  k *( (coord3[2] - 2*coord2[2] + coord1[2])*A -
                              (coord2[2] - coord1[2])*B +
                              (coord3[2] - coord2[2])*C );
            
            //force on i-1, f = k*(A*l - B*l2):
            force3[0] +=  k *( (coord2[0] - coord1[0])*A -
                              (coord3[0] - coord2[0])*C );
            
            force3[1] +=  k *( (coord2[1] - coord1[1])*A -
                              (coord3[1] - coord2[1])*C );
            
            force3[2] +=  k *( (coord2[2] - coord1[2])*A -
                              (coord3[2] - coord2[2])*C );

        }
    }
    
}



