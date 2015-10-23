
//------------------------------------------------------------------
//  **M3SYM** - Simulation Package for the Mechanochemical
//              Dynamics of Active Networks, 3rd Generation
//
//  Copyright (2014) Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the Papoian lab page for more information:
//  http://papoian.chem.umd.edu/
//------------------------------------------------------------------

#include "BubbleCylinderRepulsion.h"

#include "BubbleCylinderRepulsionExp.h"

#include "Bubble.h"
#include "Cylinder.h"
#include "Bead.h"

template <class BRepulsionInteractionType>
double BubbleCylinderRepulsion<BRepulsionInteractionType>::computeEnergy(double d) {
    
    double U = 0;
    double U_i;
    
    for (auto bb: Bubble::getBubbles()) {
        
        for(auto &c : _neighborList->getNeighbors(bb)) {
            
            double kRep = bb->getRepulsionConst();
            double screenLength = bb->getScreeningLength();
            
            double radius = bb->getRadius();
            
            Bead* bd1 = bb->getBead();
            
            //potential acts on second bead unless this is a minus end
            Bead* bd2;
            if(c->isMinusEnd())
                bd2 = c->getFirstBead();
            else
                bd2 = c->getSecondBead();
            
            if (d == 0.0)
                U_i =  _FFType.energy(bd1, bd2, radius, kRep, screenLength);
            else
                U_i = _FFType.energy(bd1, bd2, radius, kRep, screenLength, d);
            
            if(fabs(U_i) == numeric_limits<double>::infinity()
               || U_i != U_i || U_i < -1.0) {
                
                //set culprits and return
                _otherCulprit = c;
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
void BubbleCylinderRepulsion<BRepulsionInteractionType>::computeForces() {
    
    for (auto bb : Bubble::getBubbles()) {
        
        for(auto &c : _neighborList->getNeighbors(bb)) {
            
            double kRep = bb->getRepulsionConst();
            double screenLength = bb->getScreeningLength();
            
            double radius = bb->getRadius();
            
            Bead* bd1 = bb->getBead();
            
            //potential acts on second bead unless this is a minus end
            Bead* bd2;
            if(c->isMinusEnd())
                bd2 = c->getFirstBead();
            else
                bd2 = c->getSecondBead();
            
            _FFType.forces(bd1, bd2, radius, kRep, screenLength);
            
        }
    }
}


template <class BRepulsionInteractionType>
void BubbleCylinderRepulsion<BRepulsionInteractionType>::computeForcesAux() {
    
    for (auto bb : Bubble::getBubbles()) {
        
        for(auto &c : _neighborList->getNeighbors(bb)) {
            
            double kRep = bb->getRepulsionConst();
            double screenLength = bb->getScreeningLength();
            
            double radius = bb->getRadius();
            
            Bead* bd1 = bb->getBead();
            
            //potential acts on second bead unless this is a minus end
            Bead* bd2;
            if(c->isMinusEnd())
                bd2 = c->getFirstBead();
            else
                bd2 = c->getSecondBead();
            
            _FFType.forcesAux(bd1, bd2, radius, kRep, screenLength);
            
        }
    }
}

///Template specializations
template double BubbleCylinderRepulsion<BubbleCylinderRepulsionExp>::computeEnergy(double d);
template void BubbleCylinderRepulsion<BubbleCylinderRepulsionExp>::computeForces();
template void BubbleCylinderRepulsion<BubbleCylinderRepulsionExp>::computeForcesAux();
