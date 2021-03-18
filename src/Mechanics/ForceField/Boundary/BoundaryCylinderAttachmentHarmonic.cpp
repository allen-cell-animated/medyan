
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

#include "BoundaryCylinderAttachmentHarmonic.h"
#include "BoundaryCylinderAttachment.h"

#include "Bead.h"

#include "MathFunctions.h"

using namespace mathfunc;

floatingpoint BoundaryCylinderAttachmentHarmonic::energy(
    floatingpoint *coord, int *beadSet,
    floatingpoint *kattr, const std::vector< Vec< 3, floatingpoint > >& pins
) const {
    

    int n = BoundaryCylinderAttachment<BoundaryCylinderAttachmentHarmonic>::n;
    int nint = Bead::getPinnedBeads().size();

    floatingpoint U_i;
    floatingpoint U = 0;

    for(int i = 0; i < nint; i += 1) {

        const auto coord1 = makeRefVec< 3 >(coord + beadSet[n * i]);

        const auto distsq = distance2(coord1, pins[i]);
        U_i = 0.5 * kattr[i] * distsq;

        if(fabs(U_i) == numeric_limits<floatingpoint>::infinity()
           || U_i != U_i || U_i < -1.0) {

            //set culprit and return
            BoundaryInteractions::_otherCulprit = Bead::getPinnedBeads()[i];

            return -1;
        }

        U += U_i;
    }
    return U;
}

floatingpoint BoundaryCylinderAttachmentHarmonic::energy(floatingpoint *coord, floatingpoint *f, int *beadSet,
                                                  floatingpoint *kattr, floatingpoint *pins, floatingpoint d) {


    int n = BoundaryCylinderAttachment<BoundaryCylinderAttachmentHarmonic>::n;
    int nint = Bead::getPinnedBeads().size();

    floatingpoint *coord1, *pin1, *zero, dist, U_i;
    floatingpoint *force1;
    floatingpoint U = 0;
    zero = new floatingpoint[3]; zero[0] = 0; zero[1] = 0; zero[2] = 0;

    for(int i = 0; i < nint; i += 1) {

        coord1 = &coord[beadSet[n * i]];
        force1 = &f[beadSet[n * i]];

        pin1 = &pins[beadSet[n * i]];

        dist = twoPointDistanceStretched(coord1, force1, pin1, zero, d);
        U_i = 0.5 * kattr[i] * dist * dist;

        if(fabs(U_i) == numeric_limits<floatingpoint>::infinity()
           || U_i != U_i || U_i < -1.0) {

            //set culprit and return
            BoundaryInteractions::_otherCulprit = Bead::getPinnedBeads()[i];

            return -1;
        }

        U += U_i;
    }
    delete zero;
    return U;
}

void BoundaryCylinderAttachmentHarmonic::forces(
    floatingpoint *coord, floatingpoint *f, int *beadSet,
    floatingpoint *kattr, const std::vector< Vec< 3, floatingpoint > >& pins
) const {
    
    int n = BoundaryCylinderAttachment<BoundaryCylinderAttachmentHarmonic>::n;
    int nint = Bead::getPinnedBeads().size();

    for(int i = 0; i < nint; i += 1) {

        const auto coord1 = makeRefVec< 3 >(coord + beadSet[n * i]);
        auto       force1 = makeRefVec< 3 >(f     + beadSet[n * i]);

        force1 += kattr[i] * (pins[i] - coord1);
    }
}
