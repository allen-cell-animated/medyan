
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

floatingpoint MTOCAttachmentHarmonic::energy(
	floatingpoint *coord, int *beadSet, std::size_t beadStartIndex, std::size_t bubbleStartIndex,
	floatingpoint *kstr, floatingpoint* radiusvec
) const {

	floatingpoint *coord1, *coord2, dist;
	floatingpoint U = 0.0;
	floatingpoint U_i = 0.0;
	//number of beads per interaction
    int n = MTOCAttachment<MTOCAttachmentHarmonic>::n;
    //number of interactions
	int nint = MTOCAttachment<MTOCAttachmentHarmonic>::numInteractions;

    for(unsigned i = 0;i < nint; i++) {
        coord1 = &coord[beadSet[n*i]]; //coordinate of MTOC
        coord2 = &coord[beadSet[n * i + 1]];
        dist = twoPointDistance(coord1, coord2) - radiusvec[i];
        U_i = 0.5 * kstr[i] * dist * dist;

	    //set culprits and return
        if(fabs(U_i) == numeric_limits<double>::infinity()
        || U_i != U_i || U_i < -1.0) {
	        for(auto mtoc : MTOC::getMTOCs()) {
		        Bubble* b1 = mtoc->getBubble();
	        	if(b1->getIndex() * 3 + bubbleStartIndex != beadSet[n*i]) continue;
		        BubbleInteractions::_bubbleCulprit = mtoc->getBubble();
		        for (int fIndex = 0; fIndex < mtoc->getFilaments().size(); fIndex++) {
			        Filament *f = mtoc->getFilaments()[fIndex];
			        if(
						f->getMinusEndCylinder()->getFirstBead()->getIndex() * 3 + beadStartIndex
						== beadSet[n * i + 1]
					) {
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

	for(unsigned i = 0;i < nint; i++) {
		coord1 = &coord[beadSet[n*i]]; //coordinate of MTOC
		f1 = &f[beadSet[n*i]];
		coord2 = &coord[beadSet[n * i + 1]];
		f2 = &f[beadSet[n*i +1]];
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

