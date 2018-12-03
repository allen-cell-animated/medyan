
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

#include "BubbleCylinderRepulsion.h"

#include "BubbleCylinderRepulsionExp.h"

#include "MTOC.h"
#include "Bubble.h"
#include "Cylinder.h"
#include "Bead.h"

#include "MathFunctions.h"

using namespace mathfunc;

template <class BRepulsionInteractionType>
void BubbleCylinderRepulsion<BRepulsionInteractionType>::vectorize() {
    //cout << "Add later!" <<endl;
}

template <class BRepulsionInteractionType>
void BubbleCylinderRepulsion<BRepulsionInteractionType>::deallocate() {
    //cout << "Add later!" <<endl;
}

template <class BRepulsionInteractionType>
double BubbleCylinderRepulsion<BRepulsionInteractionType>::computeEnergy(double* coord, double *f, double d) {
    
    double U = 0.0;
    double U_i=0.0;
    
    for (auto bb: Bubble::getBubbles()) {
        
        for(auto &c : _neighborList->getNeighbors(bb)) {
            
            //if part of an MTOC, skip
            if(bb->isMTOC()) {
                
                auto mtoc = (MTOC*)bb->getParent();
                auto filaments = mtoc->getFilaments();
                
                auto f = (Filament*)c->getParent();
                
                if(find(filaments.begin(), filaments.end(), f) != filaments.end())
                continue;
            }
            
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
void BubbleCylinderRepulsion<BRepulsionInteractionType>::computeForces(double *coord, double *f) {
    
    for (auto bb : Bubble::getBubbles()) {
        
        for(auto &c : _neighborList->getNeighbors(bb)) {
            
            //if part of an MTOC, skip
            if(bb->isMTOC()) {
                
                auto mtoc = (MTOC*)bb->getParent();
                auto filaments = mtoc->getFilaments();
                
                auto f = (Filament*)c->getParent();
                
                if(find(filaments.begin(), filaments.end(), f) != filaments.end())
                continue;
            }
            
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


//template <class BRepulsionInteractionType>
//void BubbleCylinderRepulsion<BRepulsionInteractionType>::computeForcesAux(double *coord, double *f) {
//
//    for (auto bb : Bubble::getBubbles()) {
//
//        for(auto &c : _neighborList->getNeighbors(bb)) {
//
//            //if part of an MTOC, skip
//            if(bb->isMTOC()) {
//
//                auto mtoc = (MTOC*)bb->getParent();
//                auto filaments = mtoc->getFilaments();
//
//                auto f = (Filament*)c->getParent();
//
//                if(find(filaments.begin(), filaments.end(), f) != filaments.end())
//                    continue;
//            }
//
//            double kRep = bb->getRepulsionConst();
//            double screenLength = bb->getScreeningLength();
//
//            double radius = bb->getRadius();
//
//            Bead* bd1 = bb->getBead();
//
//            //potential acts on second bead unless this is a minus end
//            Bead* bd2;
//            if(c->isMinusEnd())
//                bd2 = c->getFirstBead();
//            else
//                bd2 = c->getSecondBead();
//
//            _FFType.forcesAux(bd1, bd2, radius, kRep, screenLength);
//
//        }
//    }
//}

template <class BRepulsionInteractionType>
void BubbleCylinderRepulsion<BRepulsionInteractionType>::computeLoadForces() {
    
    for (auto bb : Bubble::getBubbles()) {
        
        for(auto &c : _neighborList->getNeighbors(bb)) {
            
            //if part of an MTOC, skip
            if(bb->isMTOC()) {
                
                auto mtoc = (MTOC*)bb->getParent();
                auto filaments = mtoc->getFilaments();
                
                auto f = (Filament*)c->getParent();
                
                if(find(filaments.begin(), filaments.end(), f) != filaments.end())
                continue;
            }
            
            double kRep = bb->getRepulsionConst();
            double screenLength = bb->getScreeningLength();
            
            double radius = bb->getRadius();
            
            Bead* bd1 = bb->getBead();
            
            //potential acts on second bead unless this is a minus end
            Bead* bd2;
            Bead* bo;
            if(c->isPlusEnd()) {
                bd2 = c->getSecondBead();
                bo = c->getFirstBead();
                
                ///this normal is in the direction of polymerization
                auto normal = normalizeVector(twoPointDirection(bo->coordinate, bd2->coordinate));
                
                //array of coordinate values to update
                auto monSize = SysParams::Geometry().monomerSize[bd2->getType()];
                auto cylSize = SysParams::Geometry().cylinderNumMon[bd2->getType()];
                
                bd2->lfip = 0;
                for (int i = 0; i < cylSize; i++) {
                    
                    auto newCoord = vector<double>{bd2->coordinate[0] + i * normal[0] * monSize,
                        bd2->coordinate[1] + i * normal[1] * monSize,
                        bd2->coordinate[2] + i * normal[2] * monSize};
                    
                    // Projection magnitude ratio on the direction of the cylinder
                    // (Effective monomer size) = (monomer size) * proj
                    double proj = dotProduct(twoPointDirection(newCoord, bd1->coordinate), normal);
                    
                    double loadForce = _FFType.loadForces(bd1, bd2, radius, kRep, screenLength);
                    bd2->loadForcesP[bd2->lfip++] += proj * loadForce;
                }
                //reset lfi
                bd2->lfip = 0;
                
            }
            if(c->isMinusEnd()) {
                bd2 = c->getFirstBead();
                bo = c->getSecondBead();
                
                ///this normal is in the direction of polymerization
                auto normal = normalizeVector(twoPointDirection(bo->coordinate, bd2->coordinate));
                
                //array of coordinate values to update
                auto monSize = SysParams::Geometry().monomerSize[bd2->getType()];
                auto cylSize = SysParams::Geometry().cylinderNumMon[bd2->getType()];
                
                bd2->lfim = 0;
                for (int i = 0; i < cylSize; i++) {
                    
                    auto newCoord = vector<double>{bd2->coordinate[0] + i * normal[0] * monSize,
                        bd2->coordinate[1] + i * normal[1] * monSize,
                        bd2->coordinate[2] + i * normal[2] * monSize};
                    
                    // Projection magnitude ratio on the direction of the cylinder
                    // (Effective monomer size) = (monomer size) * proj
                    double proj = dotProduct(twoPointDirection(newCoord, bd1->coordinate), normal);
                    
                    double loadForce = _FFType.loadForces(bd1, bd2, radius, kRep, screenLength);
                    bd2->loadForcesM[bd2->lfim++] += proj * loadForce;
                }
                //reset lfi
                bd2->lfim = 0;
            }
            
        }
    }
}


///Template specializations
template double BubbleCylinderRepulsion<BubbleCylinderRepulsionExp>::computeEnergy(double *coord, double *f, double d);
template void BubbleCylinderRepulsion<BubbleCylinderRepulsionExp>::computeForces(double *coord, double *f);
//template void BubbleCylinderRepulsion<BubbleCylinderRepulsionExp>::computeForcesAux(double *coord, double *f);
template void BubbleCylinderRepulsion<BubbleCylinderRepulsionExp>::computeLoadForces();
template void BubbleCylinderRepulsion<BubbleCylinderRepulsionExp>::vectorize();
template void BubbleCylinderRepulsion<BubbleCylinderRepulsionExp>::deallocate();



