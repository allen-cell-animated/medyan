
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

#include "BubbleBubbleRepulsion.h"

#include "BubbleBubbleRepulsionExp.h"

#include "Bubble.h"
#include "Bead.h"

template <class BRepulsionInteractionType>
void BubbleBubbleRepulsion<BRepulsionInteractionType>::vectorize(const FFCoordinateStartingIndex& si) {
    bubbleStartIdx_ = si.bubble;
}

template <class BRepulsionInteractionType>
void BubbleBubbleRepulsion<BRepulsionInteractionType>::deallocate() {
    //cout << "Add later!" <<endl;
}


template <class BRepulsionInteractionType>
floatingpoint BubbleBubbleRepulsion<BRepulsionInteractionType>::computeEnergy(floatingpoint* coord, bool stretched) {
    
    floatingpoint U = 0.0;
    floatingpoint U_i=0.0;
    
    for (auto bb: Bubble::getBubbles()) {
        
        for(auto &bbo : _neighborList->getNeighbors(bb)) {
            
            floatingpoint kRep = bb->getRepulsionConst();
            floatingpoint screenLength = bb->getScreeningLength();
            
            floatingpoint radius1 = bb->getRadius();
            floatingpoint radius2 = bbo->getRadius();
            
            const auto bi1 = bb ->getIndex() * 3 + bubbleStartIdx_;
            const auto bi2 = bbo->getIndex() * 3 + bubbleStartIdx_;
            
            U_i = _FFType.energy(coord, bi1, bi2, radius1, radius2, kRep, screenLength);
            
            if(fabs(U_i) == numeric_limits<floatingpoint>::infinity()
               || U_i != U_i || U_i < -1.0) {
                
                //set culprits and return
                //_otherCulprit = bbo;
                //_bubbleCulprit = bb;
                
                return -1;
            }
            else
                U += U_i;
        }
    }
    
    return U;
}

template <class BRepulsionInteractionType>
void BubbleBubbleRepulsion<BRepulsionInteractionType>::computeForces(floatingpoint *coord, floatingpoint *f) {
    
    for (auto bb : Bubble::getBubbles()) {
        
        for(auto &bbo : _neighborList->getNeighbors(bb)) {
            
            floatingpoint kRep = bb->getRepulsionConst();
            floatingpoint screenLength = bb->getScreeningLength();
            
            floatingpoint radius1 = bb->getRadius();
            floatingpoint radius2 = bbo->getRadius();
            
            const auto bi1 = bb ->getIndex() * 3 + bubbleStartIdx_;
            const auto bi2 = bbo->getIndex() * 3 + bubbleStartIdx_;
            
            _FFType.forces(coord, f, bi1, bi2, radius1, radius2, kRep, screenLength);
            
        }
    }
}


///Template specializations
template floatingpoint BubbleBubbleRepulsion<BubbleBubbleRepulsionExp>::computeEnergy(floatingpoint *coord, bool stretched);
template void BubbleBubbleRepulsion<BubbleBubbleRepulsionExp>::computeForces(floatingpoint *coord, floatingpoint *f);
//template void BubbleBubbleRepulsion<BubbleBubbleRepulsionExp>::computeForcesAux(floatingpoint *coord, floatingpoint *f);
template void BubbleBubbleRepulsion<BubbleBubbleRepulsionExp>::vectorize(const FFCoordinateStartingIndex&);
template void BubbleBubbleRepulsion<BubbleBubbleRepulsionExp>::deallocate();

