
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

#include "BoundaryCylinderRepulsion.h"

#include "BoundaryCylinderRepulsionExp.h"
#include "BoundaryElement.h"

#include "Bead.h"
#include "Cylinder.h"

#include "MathFunctions.h"

using namespace mathfunc;

template <class BRepulsionInteractionType>
void BoundaryCylinderRepulsion<BRepulsionInteractionType>::vectorize() {
    
    
    //count interactions
    int nint = 0;
    for (auto be: BoundaryElement::getBoundaryElements())
        for(auto &c : _neighborList->getNeighbors(be))
            nint++;
    
    beadSet = new int[n * nint];
    krep = new double[nint];
    slen = new double[nint];
    
    
    int nbe = BoundaryElement::getBoundaryElements().size();
    int i = 0;
    int ni = 0;
    int bindex = 0;
    
    
    for (i = 0; i < nbe; i++) {

        auto be = BoundaryElement::getBoundaryElements()[i];
        auto nn = _neighborList->getNeighbors(be).size();
        
        for (ni = 0; ni < nn; ni++) {
        
            if (_neighborList->getNeighbors(be)[ni]->isPlusEnd())
                bindex = _neighborList->getNeighbors(be)[ni]->getSecondBead()->_dbIndex;
            else if(_neighborList->getNeighbors(be)[ni]->isMinusEnd())
                bindex = _neighborList->getNeighbors(be)[ni]->getFirstBead()->_dbIndex;
            
            beadSet[n * (i + ni)] = bindex;
            krep[i + ni] = be->getRepulsionConst();
            slen[i + ni] = be->getScreeningLength();
        }
    }
}

template<class BRepulsionInteractionType>
void BoundaryCylinderRepulsion<BRepulsionInteractionType>::deallocate() {
    
    delete beadSet;
    delete krep;
    delete slen;
}

template <class BRepulsionInteractionType>
double BoundaryCylinderRepulsion<BRepulsionInteractionType>::computeEnergy(double *coord, double *f, double d) {

    
    double U_i;
    
    if (d == 0.0) {
        U_i = _FFType.energy(coord, f, beadSet, krep, slen);
    }
    else {
        U_i = _FFType.energy(coord, f, beadSet, krep, slen, d);
    }
    
    return U_i;
    
    
//    double U = 0;
//    double U_i;
//    
//    for (auto be: BoundaryElement::getBoundaryElements()) {
//        
//        for(auto &c : _neighborList->getNeighbors(be)) {
//            
//            double kRep = be->getRepulsionConst();
//            double screenLength = be->getScreeningLength();
//            
//            //potential acts on second bead unless this is a minus end
//            Bead* bd;
//            if(c->isMinusEnd()) {
//                
//                bd = c->getFirstBead();
//                
//                if (d == 0.0)
//                    U_i =  _FFType.energy(bd, be->distance(bd->coordinate), kRep, screenLength);
//                else
//                    U_i = _FFType.energy(bd, be->stretchedDistance(bd->coordinate, bd->force, d), kRep, screenLength);
//                
//                if(fabs(U_i) == numeric_limits<double>::infinity()
//                   || U_i != U_i || U_i < -1.0) {
//                    
//                    //set culprits and return
//                    _otherCulprit = c;
//                    _boundaryElementCulprit = be;
//                    
//                    return -1;
//                }
//                else {
//                    U += U_i;
//                }
//                
//                
//                bd = c->getSecondBead();
//                    
//                if (d == 0.0)
//                    U_i =  _FFType.energy(bd, be->distance(bd->coordinate), kRep, screenLength);
//                else
//                    U_i = _FFType.energy(bd, be->stretchedDistance(bd->coordinate, bd->force, d), kRep, screenLength);
//                    
//                if(fabs(U_i) == numeric_limits<double>::infinity()
//                    || U_i != U_i || U_i < -1.0) {
//                        
//                    //set culprits and return
//                    _otherCulprit = c;
//                    _boundaryElementCulprit = be;
//                        
//                    return -1;
//                }
//                else {
//                    U += U_i;
//                }
//                    
//            }
//            
//            else {
//                bd = c->getSecondBead();
//                
//                if (d == 0.0)
//                    U_i =  _FFType.energy(bd, be->distance(bd->coordinate), kRep, screenLength);
//                else
//                    U_i = _FFType.energy(bd, be->stretchedDistance(bd->coordinate, bd->force, d), kRep, screenLength);
//                
//                if(fabs(U_i) == numeric_limits<double>::infinity()
//                   || U_i != U_i || U_i < -1.0) {
//                    
//                    //set culprits and return
//                    _otherCulprit = c;
//                    _boundaryElementCulprit = be;
//                    
//                    return -1;
//                }
//                else {
//                    U += U_i;
//                }
//            }
//        }
//    }
//
//    return U;
}

template <class BRepulsionInteractionType>
void BoundaryCylinderRepulsion<BRepulsionInteractionType>::computeForces(double *coord, double *f) {
    
    _FFType.forces(coord, f, beadSet, krep, slen);
    
    //    for (auto be: BoundaryElement::getBoundaryElements()) {
//        
//        for(auto &c: _neighborList->getNeighbors(be)) {
//            
//            double kRep = be->getRepulsionConst();
//            double screenLength = be->getScreeningLength();
//            
//            //potential acts on second cylinder bead unless this is a minus end
//            Bead* bd;
//            if(c->isMinusEnd()) {
//                bd = c->getFirstBead();
//                auto normal = be->normal(bd->coordinate);
//                _FFType.forces(bd, be->distance(bd->coordinate), normal, kRep, screenLength);
//                
//                bd = c->getSecondBead();
//                normal = be->normal(bd->coordinate);
//                _FFType.forces(bd, be->distance(bd->coordinate), normal, kRep, screenLength);
//            }
//            else {
//                bd = c->getSecondBead();
//                auto normal = be->normal(bd->coordinate);
//                _FFType.forces(bd, be->distance(bd->coordinate), normal, kRep, screenLength);
//            }
//            
//        }
//    }
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
                auto normal = normalizeVector(twoPointDirection(bo->coordinate, bd->coordinate));
                
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
                auto normal = normalizeVector(twoPointDirection(bo->coordinate, bd->coordinate));
                
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
template double BoundaryCylinderRepulsion<BoundaryCylinderRepulsionExp>::computeEnergy(double *coord, double *f, double d);
template void BoundaryCylinderRepulsion<BoundaryCylinderRepulsionExp>::computeForces(double *coord, double *f);
template void BoundaryCylinderRepulsion<BoundaryCylinderRepulsionExp>::computeLoadForces();

