
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

#include "BoundaryCylinderRepulsionIn.h"

#include "BoundaryCylinderRepulsionExpIn.h"
#include "BoundaryElement.h"

#include "Bead.h"
#include "Cylinder.h"

#include "MathFunctions.h"

using namespace mathfunc;

template <class BRepulsionInteractionType>
double BoundaryCylinderRepulsionIn<BRepulsionInteractionType>::computeEnergy(bool stretched) {
    
    double U = 0.0;
    double U_i=0.0;
    
    for (auto be: BoundaryElement::getBoundaryElements()) {
        
        for(auto &c : _neighborList->getNeighbors(be)) {
            
            double kRep = be->getRepulsionConst();
            double screenLength = be->getScreeningLength();
            
            //potential acts on second bead unless this is a minus end
            Bead* bd;
            if(c->isMinusEnd()) {
                
                bd = c->getFirstBead();
                
                U_i =  _FFType.energy(bd, be->distance(stretched ? bd->getCoordinate<true>() : bd->getCoordinate<false>()), kRep, screenLength);
                
                if(fabs(U_i) == numeric_limits<double>::infinity()
                   || U_i != U_i || U_i < -1.0) {
                    
                    //set culprits and return
                    _otherCulprit = c;
                    _boundaryElementCulprit = be;
                    
                    return -1;
                }
                else {
                    U += U_i;
                }
                
                
                bd = c->getSecondBead();
                
                U_i =  _FFType.energy(bd, be->distance(stretched ? bd->getCoordinate<true>() : bd->getCoordinate<false>()), kRep, screenLength);
                
                if(fabs(U_i) == numeric_limits<double>::infinity()
                   || U_i != U_i || U_i < -1.0) {
                    
                    //set culprits and return
                    _otherCulprit = c;
                    _boundaryElementCulprit = be;
                    
                    return -1;
                }
                else {
                    U += U_i;
                }
                
            }
            
            else {
                bd = c->getSecondBead();
                
                U_i =  _FFType.energy(bd, be->distance(stretched ? bd->getCoordinate<true>() : bd->getCoordinate<false>()), kRep, screenLength);
                
                if(fabs(U_i) == numeric_limits<double>::infinity()
                   || U_i != U_i || U_i < -1.0) {
                    
                    //set culprits and return
                    _otherCulprit = c;
                    _boundaryElementCulprit = be;
                    
                    return -1;
                }
                else {
                    U += U_i;
                }
            }
        }
    }
    
    return U;
}

template <class BRepulsionInteractionType>
void BoundaryCylinderRepulsionIn<BRepulsionInteractionType>::computeForces() {
    
    for (auto be: BoundaryElement::getBoundaryElements()) {
        
        for(auto &c: _neighborList->getNeighbors(be)) {
            
            double kRep = be->getRepulsionConst();
            double screenLength = be->getScreeningLength();
            
            //potential acts on second cylinder bead unless this is a minus end
            Bead* bd;
            if(c->isMinusEnd()) {
                bd = c->getFirstBead();
                auto normal = be->normal(bd->coordinate);
                _FFType.forces(bd, be->distance(bd->coordinate), normal, kRep, screenLength);
                
                bd = c->getSecondBead();
                normal = be->normal(bd->coordinate);
                _FFType.forces(bd, be->distance(bd->coordinate), normal, kRep, screenLength);
            }
            else {
                bd = c->getSecondBead();
                auto normal = be->normal(bd->coordinate);
                _FFType.forces(bd, be->distance(bd->coordinate), normal, kRep, screenLength);
            }
            
        }
    }
}


template <class BRepulsionInteractionType>
void BoundaryCylinderRepulsionIn<BRepulsionInteractionType>::computeForcesAux() {
    
    for (auto be: BoundaryElement::getBoundaryElements()) {
        
        for(auto &c : _neighborList->getNeighbors(be)) {
            
            double kRep = be->getRepulsionConst();
            double screenLength = be->getScreeningLength();
            
            //potential acts on second cylinder bead unless this is a minus end, then we compute forces
            //for both first and second bead
            Bead* bd;
            if(c->isMinusEnd()) {
                bd = c->getFirstBead();
                auto normal = be->normal(bd->coordinate);
                _FFType.forcesAux(bd, be->distance(bd->coordinate), normal, kRep, screenLength);
                
                bd = c->getSecondBead();
                normal = be->normal(bd->coordinate);
                _FFType.forcesAux(bd, be->distance(bd->coordinate), normal, kRep, screenLength);
            }
            else {
                bd = c->getSecondBead();
                auto normal = be->normal(bd->coordinate);
                _FFType.forcesAux(bd, be->distance(bd->coordinate), normal, kRep, screenLength);
            }
        }
    }
}

template <class BRepulsionInteractionType>
void BoundaryCylinderRepulsionIn<BRepulsionInteractionType>::computeLoadForces() {
    
    for (auto be: BoundaryElement::getBoundaryElements()) {
        
        for(auto &c : _neighborList->getNeighbors(be)) {
            
            double kRep = be->getRepulsionConst();
            double screenLength = be->getScreeningLength();
            
            
            //potential acts on second cylinder bead unless this is a minus end
            Bead* bd;
            Bead* bo;
            if(c->isPlusEnd()) {
                
                bd = c->getSecondBead();
                bo = c->getFirstBead();
                
                ///this normal is in the direction of polymerization
                auto normal = normalizedVector(twoPointDirection(bo->coordinate, bd->coordinate));
                
                //array of coordinate values to update
                auto monSize = SysParams::Geometry().monomerSize[bd->getType()];
                auto cylSize = SysParams::Geometry().cylinderNumMon[bd->getType()];
                
                bd->lfip = 0;
                for (int i = 0; i < cylSize; i++) {
                    
                    auto newCoord = vector<double>{bd->coordinate[0] + i * normal[0] * monSize,
                        bd->coordinate[1] + i * normal[1] * monSize,
                        bd->coordinate[2] + i * normal[2] * monSize};
                    
                    double loadForce = _FFType.loadForces(be->distance(newCoord), kRep, screenLength);
                    bd->loadForcesP[bd->lfip++] += loadForce;
                }
                //reset lfi
                bd->lfip = 0;
            }
            
            if(c->isMinusEnd()) {
                
                bd = c->getFirstBead();
                bo = c->getSecondBead();
                
                ///this normal is in the direction of polymerization
                auto normal = normalizedVector(twoPointDirection(bo->coordinate, bd->coordinate));
                
                //array of coordinate values to update
                auto monSize = SysParams::Geometry().monomerSize[bd->getType()];
                auto cylSize = SysParams::Geometry().cylinderNumMon[bd->getType()];
                
                
                bd->lfim = 0;
                for (int i = 0; i < cylSize; i++) {
                    
                    auto newCoord = vector<double>{bd->coordinate[0] + i * normal[0] * monSize,
                        bd->coordinate[1] + i * normal[1] * monSize,
                        bd->coordinate[2] + i * normal[2] * monSize};
                    
                    double loadForce = _FFType.loadForces(be->distance(newCoord), kRep, screenLength);
                    bd->loadForcesM[bd->lfim++] += loadForce;
                }
                //reset lfi
                bd->lfim = 0;
            }
            
        }
    }
}


///Template specializations
template double BoundaryCylinderRepulsionIn<BoundaryCylinderRepulsionExpIn>::computeEnergy(bool stretched);
template void BoundaryCylinderRepulsionIn<BoundaryCylinderRepulsionExpIn>::computeForces();
template void BoundaryCylinderRepulsionIn<BoundaryCylinderRepulsionExpIn>::computeForcesAux();
template void BoundaryCylinderRepulsionIn<BoundaryCylinderRepulsionExpIn>::computeLoadForces();

