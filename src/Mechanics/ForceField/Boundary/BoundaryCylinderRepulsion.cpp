
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

#include "BoundaryCylinderRepulsion.h"

#include "BoundaryCylinderRepulsionExp.h"
#include "BoundaryElement.h"

#include "Bead.h"
#include "Cylinder.h"

#include "MathFunctions.h"

using namespace mathfunc;

template <class BRepulsionInteractionType>
double BoundaryCylinderRepulsion<BRepulsionInteractionType>::computeEnergy(bool stretched) {
    
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
void BoundaryCylinderRepulsion<BRepulsionInteractionType>::computeForces() {
    
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
void BoundaryCylinderRepulsion<BRepulsionInteractionType>::computeForcesAux() {
    
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
void BoundaryCylinderRepulsion<BRepulsionInteractionType>::computeLoadForces() {
    
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
                    
                    // Projection magnitude ratio on the direction of the cylinder
                    // (Effective monomer size) = (monomer size) * proj
                    double proj = -dotProduct(be->normal(newCoord), normal);
                    if(proj < 0.0) proj = 0.0;
                
                    double loadForce = _FFType.loadForces(be->distance(newCoord), kRep, screenLength);
                    // The load force stored in bead also considers effective monomer size.
                    bd->loadForcesP[bd->lfip++] += proj * loadForce;
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
                    
                    // Projection magnitude ratio on the direction of the cylinder
                    // (Effective monomer size) = (monomer size) * proj
                    double proj = -dotProduct(be->normal(newCoord), normal);
                    if(proj < 0.0) proj = 0.0;

                    double loadForce = _FFType.loadForces(be->distance(newCoord), kRep, screenLength);
                    // The load force stored in bead also considers effective monomer size.
                    bd->loadForcesM[bd->lfim++] += proj * loadForce;
                }
                //reset lfi
                bd->lfim = 0;
            }
            
        }
    }
}


///Template specializations
template double BoundaryCylinderRepulsion<BoundaryCylinderRepulsionExp>::computeEnergy(double d);
template void BoundaryCylinderRepulsion<BoundaryCylinderRepulsionExp>::computeForces();
template void BoundaryCylinderRepulsion<BoundaryCylinderRepulsionExp>::computeForcesAux();
template void BoundaryCylinderRepulsion<BoundaryCylinderRepulsionExp>::computeLoadForces();

