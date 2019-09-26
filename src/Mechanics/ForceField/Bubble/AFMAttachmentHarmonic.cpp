
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

#include "AFMAttachmentHarmonic.h"
#include "AFMAttachment.h"
#include "Bead.h"
#include "Bubble.h"
#include "AFM.h"

#include "MathFunctions.h"

using namespace mathfunc;

double AFMAttachmentHarmonic::energy(double *coord, double *f, int *beadSet,
                                      double *kstr, double radius){
    
    double *coord1, *coord2, dist;
    double U = 0.0;
    
    int n = AFMAttachment<AFMAttachmentHarmonic>::n;
    for(auto afm : AFM::getAFMs()) {
        int nint = afm->getFilaments().size();
        
        coord1 = &coord[3 * beadSet[0]]; //coordinate of AFM
        
        for(int i = 1; i < nint + 1; i+=1){
            coord2 = &coord[3 * beadSet[n * i]];
            
            dist = twoPointDistance(coord1, coord2) - radius;
            
            U += 0.5 * kstr[i] * dist * dist;
        }
    }
    return U;
    
}

double AFMAttachmentHarmonic::energy(double *coord, double *f, int *beadSet,
                                      double *kstr, double radius, double d){
    
    double *coord1, *coord2, *f1, *f2, dist;
    double U = 0.0;
    
    int n = AFMAttachment<AFMAttachmentHarmonic>::n;
    for(auto afm : AFM::getAFMs()) {
        int nint = afm->getFilaments().size();
        
        coord1 = &coord[3 * beadSet[0]]; //coordinate of AFM
        f1 = &f[3 * beadSet[0]];
        
        for(int i = 1; i < nint + 1; i+=1){
            coord2 = &coord[3 * beadSet[n * i]];
            f2 = &f[3 * beadSet[n * i]];
            
            dist = twoPointDistanceStretched(coord1, f1,  coord2, f2, d) - radius;
            
            U += 0.5 * kstr[i] * dist * dist;
        }
    }
    return U;
    
}

void AFMAttachmentHarmonic::forces(double *coord, double *f, int *beadSet,
                                    double *kstr, double radius){
    
    
    int n = AFMAttachment<AFMAttachmentHarmonic>::n;
    for(auto afm : AFM::getAFMs()) {
        int nint = afm->getFilaments().size();
        
        double *coord1, *coord2, dist, invL;
        double f0, *f1, *f2;
        
        coord1 = &coord[3 * beadSet[0]]; //coordinate of AFM
        f1 = &f[3 * beadSet[0]];
        
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
            //            f1[0] +=  f0 * ( coord2[0] - coord1[0] );
            //            f1[1] +=  f0 * ( coord2[1] - coord1[1] );
            //            f1[2] +=  f0 * ( coord2[2] - coord1[2] );
            
        }
    }
    
}




