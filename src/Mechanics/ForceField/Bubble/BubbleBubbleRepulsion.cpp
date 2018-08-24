
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

#include "BubbleBubbleRepulsion.h"

#include "BubbleBubbleRepulsionExp.h"

#include "Bubble.h"
#include "Bead.h"

template <class BRepulsionInteractionType>
double BubbleBubbleRepulsion<BRepulsionInteractionType>::computeEnergy(double d) {
    
    double U = 0.0;
    double U_i=0.0;
    
    for (auto bb: Bubble::getBubbles()) {
        
        for(auto &bbo : _neighborList->getNeighbors(bb)) {
            
            double kRep = bb->getRepulsionConst();
            double screenLength = bb->getScreeningLength();
            
            double radius1 = bb->getRadius();
            double radius2 = bbo->getRadius();
            
            Bead* bd1 = bb->getBead();
            Bead* bd2 = bbo->getBead();
            
            if (d == 0.0)
                U_i =  _FFType.energy(
                bd1, bd2, radius1, radius2, kRep, screenLength);
            else
                U_i = _FFType.energy(
                bd1, bd2, radius1, radius2, kRep, screenLength, d);
            
            if(fabs(U_i) == numeric_limits<double>::infinity()
               || U_i != U_i || U_i < -1.0) {
                
                //set culprits and return
                _otherCulprit = bbo;
                _bubbleCulprit = bb;
                
                return -1;
            }
            else
                U += U_i;
        }
    }
    
    return U;
}

template <class BRepulsionInteractionType>
void BubbleBubbleRepulsion<BRepulsionInteractionType>::computeForces() {
    
    for (auto bb : Bubble::getBubbles()) {
        
        for(auto &bbo : _neighborList->getNeighbors(bb)) {
            
            double kRep = bb->getRepulsionConst();
            double screenLength = bb->getScreeningLength();
            
            double radius1 = bb->getRadius();
            double radius2 = bbo->getRadius();
            
            Bead* bd1 = bb->getBead();
            Bead* bd2 = bbo->getBead();
            
            _FFType.forces(bd1, bd2, radius1, radius2, kRep, screenLength);
            
        }
    }
}


template <class BRepulsionInteractionType>
void BubbleBubbleRepulsion<BRepulsionInteractionType>::computeForcesAux() {
    
    for (auto bb : Bubble::getBubbles()) {
        
        for(auto &bbo : _neighborList->getNeighbors(bb)) {
            
            double kRep = bb->getRepulsionConst();
            double screenLength = bb->getScreeningLength();
            
            double radius1 = bb->getRadius();
            double radius2 = bbo->getRadius();
            
            Bead* bd1 = bb->getBead();
            Bead* bd2 = bbo->getBead();
            
            _FFType.forcesAux(bd1, bd2, radius1, radius2, kRep, screenLength);
            
        }
    }
}

///Template specializations
template double BubbleBubbleRepulsion<BubbleBubbleRepulsionExp>::computeEnergy(double d);
template void BubbleBubbleRepulsion<BubbleBubbleRepulsionExp>::computeForces();
template void BubbleBubbleRepulsion<BubbleBubbleRepulsionExp>::computeForcesAux();
