
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

#include "BoundaryBubbleRepulsion.h"

#include "BoundaryBubbleRepulsionExp.h"
#include "BoundaryElement.h"

#include "Bubble.h"
#include "Bead.h"

template <class BRepulsionInteractionType>
floatingpoint BoundaryBubbleRepulsion<BRepulsionInteractionType>::computeEnergy(floatingpoint d) {

    totalenergyfloatingpoint U = 0.0;
    floatingpoint U_i=0.0;
    
    for (auto be: BoundaryElement::getBoundaryElements()) {
        
        for(auto &bb : _neighborList->getNeighbors(be)) {
            
            floatingpoint kRep = be->getRepulsionConst();
            floatingpoint screenLength = be->getScreeningLength();
            floatingpoint radius = bb->getRadius();
            
            Bead* bd = bb->getBead();
            
            if (d == 0.0)
                U_i =  _FFType.energy(
                bd, be->distance(bd->coordinate), radius, kRep, screenLength);
            else
                U_i = _FFType.energy(
                bd, be->stretchedDistance(bd->coordinate, bd->force, d),
                                          radius, kRep, screenLength);
            
            if(fabs(U_i) == numeric_limits<floatingpoint>::infinity()
               || U_i != U_i || U_i < -1.0) {
                
                //set culprits and return
                _otherCulprit = bb;
                _boundaryElementCulprit = be;
                
                return -1;
            }
            else
                U += U_i;
        }
    }
    
    return U;
}

template <class BRepulsionInteractionType>
void BoundaryBubbleRepulsion<BRepulsionInteractionType>::computeForces() {
    
    for (auto be: BoundaryElement::getBoundaryElements()) {
        
        for(auto &bb : _neighborList->getNeighbors(be)) {
            
            floatingpoint kRep = be->getRepulsionConst();
            floatingpoint screenLength = be->getScreeningLength();
            floatingpoint radius = bb->getRadius();
            
            Bead* bd = bb->getBead();
            
            auto normal = be->normal(bd->coordinate);
            _FFType.forces(bd, be->distance(bd->coordinate),
                           radius, normal, kRep, screenLength);
            
        }
    }
}


template <class BRepulsionInteractionType>
void BoundaryBubbleRepulsion<BRepulsionInteractionType>::computeForcesAux() {
    
    for (auto be: BoundaryElement::getBoundaryElements()) {
        
        for(auto &bb : _neighborList->getNeighbors(be)) {
            
            floatingpoint kRep = be->getRepulsionConst();
            floatingpoint screenLength = be->getScreeningLength();
            floatingpoint radius = bb->getRadius();
            
            Bead* bd = bb->getBead();
            
            auto normal = be->normal(bd->coordinate);
            _FFType.forcesAux(bd, be->distance(bd->coordinate),
                              radius, normal, kRep, screenLength);
            
        }
    }
}

///Template specializations
template floatingpoint BoundaryBubbleRepulsion<BoundaryBubbleRepulsionExp>::computeEnergy(floatingpoint d);
template void BoundaryBubbleRepulsion<BoundaryBubbleRepulsionExp>::computeForces();
template void BoundaryBubbleRepulsion<BoundaryBubbleRepulsionExp>::computeForcesAux();
