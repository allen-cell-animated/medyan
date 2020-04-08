
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

floatingpoint AFMAttachmentHarmonic::energy(
    floatingpoint *coord,
    int numInteractions, int *beadSet, floatingpoint *kstr, const floatingpoint* radii
) const {
    
    floatingpoint *coord1, *coord2, dist;
    floatingpoint U = 0.0;

    for(int i = 0; i < numInteractions; ++i) {
        coord1 = &coord[beadSet[2*i    ]]; //coordinate of AFM
        coord2 = &coord[beadSet[2*i + 1]]; //coordinate of the filament

        dist = twoPointDistance(coord1, coord2) - radii[i];

        U += 0.5 * kstr[i] * dist * dist;

    }
    return U;
    
}

void AFMAttachmentHarmonic::forces(
    floatingpoint *coord, floatingpoint* f,
    int numInteractions, int *beadSet, floatingpoint *kstr, const floatingpoint* radii
) const {
    for(int i = 0; i < numInteractions; ++i) {

        floatingpoint *coord1, *coord2, dist, invL;
        floatingpoint f0, *f1, *f2;
        
        coord1 = &coord[beadSet[2*i    ]]; //coordinate of AFM
        coord2 = &coord[beadSet[2*i + 1]]; //coordinate of the filament
        f1 = &f[beadSet[2*i    ]];
        f2 = &f[beadSet[2*i + 1]];
    
        dist = twoPointDistance(coord1, coord2);
        invL = 1 / dist;
        
        f0 = kstr[i] * ( dist - radii[i] ) * invL;
        
        // f1[0] +=  f0 * ( coord2[0] - coord1[0] );
        // f1[1] +=  f0 * ( coord2[1] - coord1[1] );
        // f1[2] +=  f0 * ( coord2[2] - coord1[2] );

        f2[0] +=  f0 * ( coord1[0] - coord2[0] );
        f2[1] +=  f0 * ( coord1[1] - coord2[1] );
        f2[2] +=  f0 * ( coord1[2] - coord2[2] );
        
        
    }
    
}




