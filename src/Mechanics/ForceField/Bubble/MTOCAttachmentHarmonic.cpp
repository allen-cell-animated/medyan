
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

#include "MTOCAttachmentHarmonic.h"
#include "MTOCAttachment.h"
#include "Bead.h"
#include "Bubble.h"
#include "MTOC.h"

#include "MathFunctions.h"

using namespace mathfunc;

double MTOCAttachmentHarmonic::energy(double *coord, double *f, int *beadSet,
                                      double *kstr, double radius){
    
    double *coord1, *coord2, dist;
    double U = 0.0;
    
    int n = MTOCAttachment<MTOCAttachmentHarmonic>::n;
    for(auto mtoc : MTOC::getMTOCs()) {
        int nint = mtoc->getFilaments().size();
        
        coord1 = &coord[0]; //coordinate of MTOC
        
        for(int i = 1; i < nint + 1; i+=1){
            coord2 = &coord[3 * beadSet[n * i]];
            
            dist = twoPointDistance(coord1, coord2) - radius;
            
            U += 0.5 * kstr[i] * dist * dist;
        }
    }
    return U;
    
}

double MTOCAttachmentHarmonic::energy(double *coord, double *f, int *beadSet,
                                      double *kstr, double radius, double d){
    
    double *coord1, *coord2, *f1, *f2, dist;
    double U = 0.0;
    
    int n = MTOCAttachment<MTOCAttachmentHarmonic>::n;
    for(auto mtoc : MTOC::getMTOCs()) {
        int nint = mtoc->getFilaments().size();
        
        coord1 = &coord[0]; //coordinate of MTOC
        f1 = &f[0];
        
        for(int i = 1; i < nint + 1; i+=1){
            coord2 = &coord[3 * beadSet[n * i]];
            f2 = &f[3 * beadSet[n * i]];
            
            dist = twoPointDistanceStretched(coord1, f1,  coord2, f2, d) - radius;
            
            U += 0.5 * kstr[i] * dist * dist;
        }
    }
    return U;
    
}

void MTOCAttachmentHarmonic::forces(double *coord, double *f, int *beadSet,
                                    double *kstr, double radius){
    
    
    int n = MTOCAttachment<MTOCAttachmentHarmonic>::n;
    for(auto mtoc : MTOC::getMTOCs()) {
        int nint = mtoc->getFilaments().size();
        
        double *coord1, *coord2, dist, invL;
        double f0, *f1, *f2;
        
        coord1 = &coord[0]; //coordinate of MTOC
        f1 = &f[0];
        
        for(int i = 1; i < nint + 1; i+=1){
            coord2 = &coord[3 * beadSet[n * i]];
            f2 = &f[3 * beadSet[n * i]];
            
            dist = twoPointDistance(coord1, coord2);
            invL = 1 / dist;
            
            f0 = kstr[i] * ( dist - radius ) * invL;
            
            f2[0] +=  f0 * ( coord1[0] - coord2[0] );
            f2[1] +=  f0 * ( coord1[1] - coord2[1] );
            f2[2] +=  f0 * ( coord1[2] - coord2[2] );
            
            // force i-1
            f1[0] +=  f0 * ( coord2[0] - coord1[0] );
            f1[1] +=  f0 * ( coord2[1] - coord1[1] );
            f1[2] +=  f0 * ( coord2[2] - coord1[2] );
            
        }
    }
    
}

//void MTOCAttachmentHarmonic::forcesAux(Bead* b1, Bead* b2, double kStretch, double radius){
//
//    cout << "MTOCAttachmentHarmonic::forcesAux should not be used in vectorized version." << endl;
//
//    double dist = twoPointDistance(b1->coordinate, b2->coordinate);
//    double invL = 1 / dist;
//
//    double f0 = kStretch * ( dist - radius ) * invL;
//
//    //force on i
//    b2->forceAux[0] +=  f0 * ( b1->coordinate[0] - b2->coordinate[0] );
//    b2->forceAux[1] +=  f0 * ( b1->coordinate[1] - b2->coordinate[1] );
//    b2->forceAux[2] +=  f0 * ( b1->coordinate[2] - b2->coordinate[2] );
//
//    // force i-1
//    b1->forceAux[0] +=  f0 * ( b2->coordinate[0] - b1->coordinate[0] );
//    b1->forceAux[1] +=  f0 * ( b2->coordinate[1] - b1->coordinate[1] );
//    b1->forceAux[2] +=  f0 * ( b2->coordinate[2] - b1->coordinate[2] );
//}


