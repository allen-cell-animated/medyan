
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v4.0
//
//  Copyright (2015-2018)  Papoian Lab, University of Maryland
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
#include "Filament.h"
#include "Cylinder.h"
#include "MathFunctions.h"

using namespace mathfunc;

floatingpoint MTOCAttachmentHarmonic::energy(floatingpoint *coord, floatingpoint *f, int *beadSet,
                                             floatingpoint *kstr, floatingpoint* radiusvec){

	floatingpoint *coord1, *coord2, dist;
	floatingpoint U = 0.0;
	floatingpoint U_i = 0.0;
	//number of beads per interaction
    int n = MTOCAttachment<MTOCAttachmentHarmonic>::n;
    //number of interactions
	int nint = MTOCAttachment<MTOCAttachmentHarmonic>::numInteractions;

    for(uint i = 0;i < nint; i++) {
        coord1 = &coord[3 * beadSet[n*i]]; //coordinate of MTOC
        coord2 = &coord[3 * beadSet[n * i + 1]];
        dist = twoPointDistance(coord1, coord2) - radiusvec[i];
        U_i = 0.5 * kstr[i] * dist * dist;

	    //set culprits and return
        if(fabs(U_i) == numeric_limits<double>::infinity()
        || U_i != U_i || U_i < -1.0) {
	        for(auto mtoc : MTOC::getMTOCs()) {
		        Bead* b1 = mtoc->getBubble()->getBead();
	        	if(b1->getStableIndex() != beadSet[n*i]) continue;
		        BubbleInteractions::_bubbleCulprit = mtoc->getBubble();
		        for (int fIndex = 0; fIndex < mtoc->getFilaments().size(); fIndex++) {
			        Filament *f = mtoc->getFilaments()[fIndex];
			        if(f->getMinusEndCylinder()->getFirstBead()->getStableIndex() == beadSet[n *
			        i + 1]){
				        BubbleInteractions::_otherCulprit = f;
				        break;
			        }
		        }
	        }
	        return -1;
        } else
        	U += U_i;
    }
    return U;
}

floatingpoint MTOCAttachmentHarmonic::energy(floatingpoint *coord, floatingpoint *f, int *beadSet,
                                             floatingpoint *kstr, floatingpoint* radiusvec,
                                             floatingpoint d){

	floatingpoint *coord1, *coord2, *f1, *f2, dist;
	floatingpoint U = 0.0;
	floatingpoint U_i = 0.0;
	//number of beads per interaction
	int n = MTOCAttachment<MTOCAttachmentHarmonic>::n;
	//number of interactions
	int nint = MTOCAttachment<MTOCAttachmentHarmonic>::numInteractions;

	for(uint i = 0;i < nint; i++) {
		coord1 = &coord[3 * beadSet[n*i]]; //coordinate of MTOC
		f1 = &f[3 * beadSet[n*i]];
		coord2 = &coord[3 * beadSet[n * i + 1]];
		f2 = &f[3 * beadSet[n*i +1]];
		dist = twoPointDistanceStretched(coord1, f1,  coord2, f2, d) - radiusvec[i];
		U_i = 0.5 * kstr[i] * dist * dist;

		//set culprits and return
		if(fabs(U_i) == numeric_limits<double>::infinity()
		   || U_i != U_i || U_i < -1.0) {
			for(auto mtoc : MTOC::getMTOCs()) {
				Bead* b1 = mtoc->getBubble()->getBead();
				if(b1->getStableIndex() != beadSet[n*i]) continue;
				BubbleInteractions::_bubbleCulprit = mtoc->getBubble();
				for (int fIndex = 0; fIndex < mtoc->getFilaments().size(); fIndex++) {
					Filament *f = mtoc->getFilaments()[fIndex];
					if(f->getMinusEndCylinder()->getFirstBead()->getStableIndex() == beadSet[n * i + 1]){
						BubbleInteractions::_otherCulprit = f;
						break;
					}
				}
			}
			return -1;
		} else
			U += U_i;
	}
    return U;

}

void MTOCAttachmentHarmonic::forces(floatingpoint *coord, floatingpoint *f, int *beadSet,
                                    floatingpoint *kstr, floatingpoint* radiusvec){

	//number of beads per interaction
	int n = MTOCAttachment<MTOCAttachmentHarmonic>::n;
	//number of interactions
	int nint = MTOCAttachment<MTOCAttachmentHarmonic>::numInteractions;
	floatingpoint *coord1, *coord2, dist, invL;
	floatingpoint f0, *f1, *f2;

	for(uint i = 0;i < nint; i++) {
		coord1 = &coord[3 * beadSet[n*i]]; //coordinate of MTOC
		f1 = &f[3 * beadSet[n*i]];
		coord2 = &coord[3 * beadSet[n * i + 1]];
		f2 = &f[3 * beadSet[n*i +1]];
		dist = twoPointDistance(coord1, coord2);
		invL = 1 / dist;

		f0 = kstr[i] * ( dist - radiusvec[i] ) * invL;

		f2[0] +=  f0 * ( coord1[0] - coord2[0] );
		f2[1] +=  f0 * ( coord1[1] - coord2[1] );
		f2[2] +=  f0 * ( coord1[2] - coord2[2] );

		// force i-1
		f1[0] +=  f0 * ( coord2[0] - coord1[0] );
		f1[1] +=  f0 * ( coord2[1] - coord1[1] );
		f1[2] +=  f0 * ( coord2[2] - coord1[2] );
	}
}

