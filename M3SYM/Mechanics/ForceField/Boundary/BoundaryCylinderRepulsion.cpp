
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

#include "BoundaryCylinderRepulsion.h"

#include "BoundaryCylinderRepulsionExp.h"
#include "BoundaryElement.h"

#include "Bead.h"
#include "Cylinder.h"

template <class BRepulsionInteractionType>
double BoundaryCylinderRepulsion<BRepulsionInteractionType>::computeEnergy(double d) {
    
    double U = 0;
    double U_i;
    
    for (auto be: BoundaryElement::getBoundaryElements()) {
        
        for(auto &c : _neighborList->getNeighbors(be)) {
            
            double kRep = be->getRepulsionConst();
            double screenLength = be->getScreeningLength();
            
            //potential acts on second bead unless this is a minus end
            Bead* bd;
            if(c->isMinusEnd())
                bd = c->getFirstBead();
            else
                bd = c->getSecondBead();
            
            if (d == 0.0)
                U_i =  _FFType.energy(
                bd, be->distance(bd->coordinate), kRep, screenLength);
            else
                U_i = _FFType.energy(
                bd, be->stretchedDistance(bd->coordinate, bd->force, d), kRep, screenLength);
            
            if(fabs(U_i) == numeric_limits<double>::infinity()
               || U_i != U_i || U_i < -1.0) {
                
                //set culprits and return
                _otherCulprit = c;
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
void BoundaryCylinderRepulsion<BRepulsionInteractionType>::computeForces() {
    
    for (auto be: BoundaryElement::getBoundaryElements()) {
        
        for(auto &c: _neighborList->getNeighbors(be)) {
            
            double kRep = be->getRepulsionConst();
            double screenLength = be->getScreeningLength();
            
            //potential acts on second cylinder bead unless this is a minus end
            Bead* bd;
            if(c->isMinusEnd())
                bd = c->getFirstBead();
            else
                bd = c->getSecondBead();
            
            auto normal = be->normal(bd->coordinate);
            _FFType.forces(bd, be->distance(bd->coordinate), normal, kRep, screenLength);
            
        }
    }
}


template <class BRepulsionInteractionType>
void BoundaryCylinderRepulsion<BRepulsionInteractionType>::computeForcesAux() {
    
    for (auto be: BoundaryElement::getBoundaryElements()) {
        
        for(auto &c : _neighborList->getNeighbors(be)) {
            
            double kRep = be->getRepulsionConst();
            double screenLength = be->getScreeningLength();
            
            //potential acts on second cylinder bead unless this is a minus end
            Bead* bd;
            if(c->isMinusEnd())
                bd = c->getFirstBead();
            else
                bd = c->getSecondBead();
            
            auto normal = be->normal(bd->coordinate);
            _FFType.forcesAux(bd, be->distance(bd->coordinate), normal, kRep, screenLength);
            
        }
    }
}

///Template specializations
template double BoundaryCylinderRepulsion<BoundaryCylinderRepulsionExp>::computeEnergy(double d);
template void BoundaryCylinderRepulsion<BoundaryCylinderRepulsionExp>::computeForces();
template void BoundaryCylinderRepulsion<BoundaryCylinderRepulsionExp>::computeForcesAux();

