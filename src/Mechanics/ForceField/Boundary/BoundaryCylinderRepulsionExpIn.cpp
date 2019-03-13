
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

#include "BoundaryCylinderRepulsionExpIn.h"
#include "BoundaryCylinderRepulsion.h"

#include "BoundaryElement.h"
#include "Bead.h"

floatingpoint BoundaryCylinderRepulsionExpIn::energy(Bead* b, floatingpoint r,
                                            floatingpoint kRep, floatingpoint screenLength) {
    //Boundary repulsion occurs 100 nm below the actual boundary
    floatingpoint R = -r/screenLength + 100/screenLength;
    return kRep * exp(R);
}

void BoundaryCylinderRepulsionExpIn::forces(Bead* b, floatingpoint r, vector<floatingpoint>& norm,
                                          floatingpoint kRep, floatingpoint screenLength) {
    
    floatingpoint R = -r/screenLength + 100/screenLength;
    floatingpoint f0 = kRep * exp(R)/screenLength;
    
    b->force[0] += f0 *norm[0];
    b->force[1] += f0 *norm[1];
    b->force[2] += f0 *norm[2];
    // Qin add brforce
    b->brforce[0] = f0 *norm[0];
    b->brforce[1] = f0 *norm[1];
    b->brforce[2] = f0 *norm[2];
    
}

void BoundaryCylinderRepulsionExpIn::forcesAux(Bead* b, floatingpoint r, vector<floatingpoint>& norm,
                                             floatingpoint kRep, floatingpoint screenLength) {
    
    floatingpoint R = -r/screenLength + 100/screenLength;
    floatingpoint f0 = kRep * exp(R)/screenLength;
    
    b->forceAux[0] += f0 *norm[0];
    b->forceAux[1] += f0 *norm[1];
    b->forceAux[2] += f0 *norm[2];
    
}

floatingpoint BoundaryCylinderRepulsionExpIn::loadForces(floatingpoint r, floatingpoint kRep, floatingpoint screenLength) {
    
    floatingpoint R = -r/screenLength + 100/screenLength;
    return kRep * exp(R)/screenLength;
    
}

floatingpoint BoundaryCylinderRepulsionExpIn::energy(floatingpoint *coord, floatingpoint *f, int *beadSet,
                                            floatingpoint *krep, floatingpoint *slen, int *nneighbors) {

    int nb, nc;
    floatingpoint *coord1, R, r, U_i;
    floatingpoint U = 0.0;
    int Cumnc=0;
    auto beList = BoundaryElement::getBoundaryElements();
    nb = beList.size();

    for (int ib = 0; ib < nb; ib++) {

        auto be = beList[ib];
        nc = nneighbors[ib];

        for (int ic = 0; ic < nc; ic++) {

            coord1 = &coord[3 * beadSet[Cumnc + ic]];
            r = be->distance(coord1);

            R = -r / slen[Cumnc + ic] + 100/slen[Cumnc + ic];
            U_i = krep[Cumnc + ic] * exp(R);

            if (fabs(U_i) == numeric_limits<floatingpoint>::infinity()
                || U_i != U_i || U_i < -1.0) {

                //set culprit and return
                BoundaryInteractions::_boundaryElementCulprit = be;
                ///TODO
                //BoundaryInteractions::_otherCulprit;

                return -1;
            }
            U += U_i;
        }
        Cumnc += nc;
    }
    return U;
}

floatingpoint BoundaryCylinderRepulsionExpIn::energy(floatingpoint *coord, floatingpoint *f, int *beadSet,
                                            floatingpoint *krep, floatingpoint *slen, int *nneighbors, floatingpoint d) {

    int nb, nc;
    floatingpoint *coord1, *force1, R, r, U_i;
    floatingpoint U = 0.0;
    int Cumnc=0;
    auto beList = BoundaryElement::getBoundaryElements();
    nb = beList.size();

    for (int ib = 0; ib < nb; ib++) {

        auto be = beList[ib];
        nc = nneighbors[ib];

        for(int ic = 0; ic < nc; ic++) {

            coord1 = &coord[3 * beadSet[Cumnc + ic]];
            force1 = &f[3 * beadSet[Cumnc + ic]];

            r = be->stretchedDistance(coord1, force1, d);

            R = -r / slen[Cumnc + ic] + 100/slen[Cumnc + ic];

            U_i = krep[Cumnc + ic] * exp(R);

            if(fabs(U_i) == numeric_limits<floatingpoint>::infinity()
               || U_i != U_i || U_i < -1.0) {

                //set culprit and return
                BoundaryInteractions::_boundaryElementCulprit = be;
                ///TODO
                //BoundaryInteractions::_otherCulprit;

                return -1;
            }
            U += U_i;
        }
        Cumnc+=nc;
    }
    return U;
}



void BoundaryCylinderRepulsionExpIn::forces(floatingpoint *coord, floatingpoint *f, int *beadSet,
                                          floatingpoint *krep, floatingpoint *slen, int *nneighbors) {
    int nb, nc;
    floatingpoint *coord1, *force1, R, r, f0;
    floatingpoint *F_i;

    auto beList = BoundaryElement::getBoundaryElements();
    nb = beList.size();
    int Cumnc=0;

    for (int ib = 0; ib < nb; ib++) {

        auto be = beList[ib];
        nc = nneighbors[ib];
        for(int ic = 0; ic < nc; ic++) {
            coord1 = &coord[3 * beadSet[ Cumnc + ic]];
            force1 = &f[3 * beadSet[ Cumnc + ic]];
            r = be->distance(coord1);
            auto norm = be->normal(coord1);

            R = -r / slen[Cumnc + ic] + 100/slen[Cumnc + ic];
            f0 = krep[Cumnc + ic] * exp(R)/ slen[Cumnc + ic];
            force1[0] += f0 *norm[0];
            force1[1] += f0 *norm[1];
            force1[2] += f0 *norm[2];

        }
        Cumnc+=nc;
    }}
