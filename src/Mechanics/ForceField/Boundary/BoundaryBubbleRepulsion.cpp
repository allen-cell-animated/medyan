
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

#include "BoundaryBubbleRepulsion.h"

#include "BoundaryBubbleRepulsionExp.h"
#include "BoundaryElement.h"

#include "Bubble.h"
#include "Bead.h"

template <class BRepulsionInteractionType>
void BoundaryBubbleRepulsion<BRepulsionInteractionType>::vectorize(const FFCoordinateStartingIndex& si) {
    //count interactions
    nint = 0;
    for (auto be: BoundaryElement::getBoundaryElements())
    {

        nint += _neighborList->getNeighbors(be).size();
    }

    beadSet = new int[n * nint];
    krep = new floatingpoint[nint];
    slen = new floatingpoint[nint];
    auto beList = BoundaryElement::getBoundaryElements();

    int nbe = BoundaryElement::getBoundaryElements().size();
    int i = 0;
    int ni = 0;
    int bindex = 0;

    nneighbors = new int[nbe];//stores number of interactions per boundary element.

    int *nintvec;

    nintvec = new int[nbe];//stores cumulative number of nneighbors.

    int cumnn=0;
    for (i = 0; i < nbe; i++) {

        auto be = BoundaryElement::getBoundaryElements()[i];//beList[i];
        auto nn = _neighborList->getNeighbors(be).size();

        nneighbors[i] = 0;
        auto idx = 0;

        for (ni = 0; ni < nn; ni++) {

            bindex = _neighborList->getNeighbors(be)[ni]->getIndex() * 3 + si.bubble;
            beadSet[cumnn+idx] = bindex;
            krep[cumnn+idx] = be->getRepulsionConst();
            slen[cumnn+idx] = be->getScreeningLength();
            idx++;


        }
        nneighbors[i]=idx;
        cumnn+=idx;
        nintvec[i] = cumnn;
    }


}

template <class BRepulsionInteractionType>
void BoundaryBubbleRepulsion<BRepulsionInteractionType>::deallocate() {
    delete [] beadSet;
    delete [] krep;
    delete [] slen;
    delete [] nneighbors;
}

template <class BRepulsionInteractionType>
floatingpoint BoundaryBubbleRepulsion<BRepulsionInteractionType>::computeEnergy(floatingpoint *coord) {
    
    return _FFType.energy(coord, beadSet, krep, slen, nneighbors);
    
}

template <class BRepulsionInteractionType>
void BoundaryBubbleRepulsion<BRepulsionInteractionType>::computeForces(floatingpoint *coord, floatingpoint *f) {
    
    _FFType.forces(coord, f, beadSet, krep, slen, nneighbors);
}


///Template specializations
template floatingpoint BoundaryBubbleRepulsion<BoundaryBubbleRepulsionExp>::computeEnergy(floatingpoint *coord);
template void BoundaryBubbleRepulsion<BoundaryBubbleRepulsionExp>::computeForces(floatingpoint *coord, floatingpoint *f);
template void BoundaryBubbleRepulsion<BoundaryBubbleRepulsionExp>::vectorize(const FFCoordinateStartingIndex&);
template void BoundaryBubbleRepulsion<BoundaryBubbleRepulsionExp>::deallocate();

