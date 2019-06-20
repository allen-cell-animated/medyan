
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v3.2.1
//
//  Copyright (2015-2018)  Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the MEDYAN web page for more information:
//  http://www.medyan.org
//------------------------------------------------------------------

#include "BubbleCylinderRepulsionExp.h"
#include "BubbleCylinderRepulsion.h"

#include "Bead.h"
#include "Bubble.h"

#include "MathFunctions.h"

using namespace mathfunc;

floatingpoint BubbleCylinderRepulsionExp::energy(floatingpoint *coord, floatingpoint *f, int *beadSet, int *bubbleSet,
                                                 floatingpoint *krep, floatingpoint *slen, floatingpoint *radius, int *nneighbors) {
    
    int nb, nc;
	floatingpoint *coord1, *coordb, R, r, U_i;
	floatingpoint U = 0.0;
    int Cumnc=0;
    auto bbList = Bubble::getBubbles();
    nb = bbList.size();
    //loop through bubbles
    for (int ib = 0; ib < nb; ib++) {

        coordb = &coord[3 * bubbleSet[ib]];
        auto be = bbList[ib];
        nc = nneighbors[ib];
        auto bradius = radius[ib];

        for (int ic = 0; ic < nc; ic++) {

            coord1 = &coord[3 * beadSet[Cumnc + ic]];
	        floatingpoint dist = twoPointDistance(coordb, coord1);
	        floatingpoint effd = dist - bradius;


            R = -effd / slen[Cumnc + ic];
            U_i = krep[Cumnc + ic] * exp(R);

            if (fabs(U_i) == numeric_limits<floatingpoint>::infinity()
                || U_i != U_i || U_i < -1.0) {

                //set culprit and return
                BubbleInteractions::_bubbleCulprit = be;

                return -1;
            }
            U += U_i;
        }
        Cumnc += nc;
    }
    return U;

}

floatingpoint BubbleCylinderRepulsionExp::energy(floatingpoint *coord, floatingpoint *f, int *beadSet, int *bubbleSet,
                                          floatingpoint *krep, floatingpoint *slen, floatingpoint *radius, int *nneighbors, floatingpoint d) {
    
    int nb, nc;
	floatingpoint *coord1, *coordb, *fb, *f1, R, r, U_i;
	floatingpoint U = 0.0;
    int Cumnc=0;
    auto bbList = Bubble::getBubbles();
    nb = bbList.size();
    //loop through bubbles
    for (int ib = 0; ib < nb; ib++) {

        coordb = &coord[3 * bubbleSet[ib]];
        fb = &f[3 * bubbleSet[ib]];
        auto be = bbList[ib];
        nc = nneighbors[ib];
        auto bradius = radius[ib];

        for (int ic = 0; ic < nc; ic++) {

            coord1 = &coord[3 * beadSet[Cumnc + ic]];
            f1 = &f[3 * beadSet[Cumnc + ic]];
//            double dist = twoPointDistanceStretched(b1->vcoordinate(), b1->force,
//                                                    b2->vcoordinate(), b2->force, d);
            floatingpoint dist = twoPointDistanceStretched(coordb, fb, coord1, f1, d);
            floatingpoint effd = dist - bradius;

            R = -effd / slen[Cumnc + ic];
            U_i = krep[Cumnc + ic] * exp(R);

            if (fabs(U_i) == numeric_limits<floatingpoint>::infinity()
                || U_i != U_i || U_i < -1.0) {

                //set culprit and return
                BubbleInteractions::_bubbleCulprit = be;

                return -1;
            }
            U += U_i;
        }
        Cumnc += nc;
    }
    return U;
}

void BubbleCylinderRepulsionExp::forces(floatingpoint *coord, floatingpoint *f, int *beadSet, int *bubbleSet,
                                        floatingpoint *krep, floatingpoint *slen, floatingpoint *radius, int *nneighbors) {

    
    //get norm
//    auto norm = normalizeVector(twoPointDirection(b1->vcoordinate(), b2->vcoordinate()));
    
    int nb, nc;
	floatingpoint *coord1, *coordb, *fb, *f1, R, f0, invL;
    int Cumnc=0;
    auto bbList = Bubble::getBubbles();
    nb = bbList.size();
    //loop through bubbles
    for (int ib = 0; ib < nb; ib++) {

        coordb = &coord[3 * bubbleSet[ib]];
        fb = &f[3 * bubbleSet[ib]];
        nc = nneighbors[ib];
        auto bradius = radius[ib];


        for (int ic = 0; ic < nc; ic++) {

            coord1 = &coord[3 * beadSet[Cumnc + ic]];
            f1 = &f[3 * beadSet[Cumnc + ic]];
	        floatingpoint dist = twoPointDistance(coordb, coord1);
            invL = 1 / dist;
	        floatingpoint effd = dist - bradius;

            R = -effd / slen[Cumnc + ic];
            f0 = krep[Cumnc + ic] * exp(R)/ slen[Cumnc + ic] * invL;

            fb[0] +=  f0 * ( coordb[0] - coord1[0] );
            fb[1] +=  f0 * ( coordb[1] - coord1[1] );
            fb[2] +=  f0 * ( coordb[2] - coord1[2] );

            f1[0] +=  f0 * ( coord1[0] - coordb[0] );
            f1[1] +=  f0 * ( coord1[1] - coordb[1] );
            f1[2] +=  f0 * ( coord1[2] - coordb[2] );

        }
        Cumnc += nc;
    }
    


}

//void BubbleCylinderRepulsionExp::forcesAux(Bead* b1, Bead* b2, double radius,
//                                           double kRep, double screenLength) {
//
//    //get dist
//    double dist = twoPointDistance(b1->coordinate, b2->coordinate);
//
//    double effd = dist - radius;
//
//    double R = -effd / screenLength;
//    double f0 = kRep * exp(R) / screenLength;
//
//    //get norm
//    auto norm = normalizeVector(twoPointDirection(b1->coordinate, b2->coordinate));
//
//    b1->force[0] += - f0 *norm[0];
//    b1->force[1] += - f0 *norm[1];
//    b1->force[2] += - f0 *norm[2];
//
//    b2->force[0] += f0 *norm[0];
//    b2->force[1] += f0 *norm[1];
//    b2->force[2] += f0 *norm[2];
//
//}

floatingpoint BubbleCylinderRepulsionExp::loadForces(Bead* b1, Bead* b2, floatingpoint radius,
                                              floatingpoint kRep, floatingpoint screenLength) {
    
    //get dist
    floatingpoint dist = twoPointDistance(b1->vcoordinate(), b2->vcoordinate());
    
    floatingpoint effd = dist - radius;
    
    floatingpoint R = -effd / screenLength;
    return kRep * exp(R) / screenLength;
}
