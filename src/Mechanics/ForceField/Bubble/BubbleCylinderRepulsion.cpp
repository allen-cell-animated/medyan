
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
    //count interactions
    nint = 0;
    for (auto bb:Bubble::getBubbles()){
        for(auto &c : _neighborList->getNeighbors(bb)) {
            //if part of an MTOC, skip
            if(bb->isMTOC()) {
                
                auto mtoc = (MTOC*)bb->getParent();
                auto filaments = mtoc->getFilaments();
                
                auto f = (Filament*)c->getParent();
                
                if(find(filaments.begin(), filaments.end(), f) != filaments.end())
                    continue;
            }
            if(c->isMinusEnd()) nint++;
            nint++;
        }

    }

    //stores number of interactions per bubble
    nneighbors = new int[Bubble::getBubbles().size()];
    //stores bubble index
    bubbleSet = new int[Bubble::getBubbles().size()];
    radius = new floatingpoint[Bubble::getBubbles().size()];
    //stores cumulative number of nneighbors, for CUDA only.
//    nintvec = new int[Bubble::getBubbles().size()];
    beadSet = new int[nint];
    krep = new floatingpoint[nint];
    slen = new floatingpoint[nint];

    int idb = 0;

    int bindex = 0;
    int cumnn=0;


    for (auto bb:Bubble::getBubbles()){
        
        nneighbors[idb] = 0;
        int idx = 0;

        //total number of neighbor cylinders
        int cmax = _neighborList->getNeighbors(bb).size();
        for(int ni = 0; ni < cmax; ni++){
            //if part of an MTOC, skip
            if(bb->isMTOC()) {
                
                auto mtoc = (MTOC*)bb->getParent();
                auto filaments = mtoc->getFilaments();
                
                auto f = (Filament*)_neighborList->getNeighbors(bb)[ni]->getParent();
                
                if(find(filaments.begin(), filaments.end(), f) != filaments.end())
                    continue;
            }
            //if this neighbor cylinder contains a minusend, add the frist bead
            if(_neighborList->getNeighbors(bb)[ni]->isMinusEnd())
            {
                bindex = _neighborList->getNeighbors(bb)[ni]->getFirstBead()->getStableIndex();
                beadSet[cumnn+idx] = bindex;
                krep[cumnn+idx] = bb->getRepulsionConst();
                slen[cumnn+idx] = bb->getScreeningLength();
                idx++;
            }
            //add all second beads
            bindex = _neighborList->getNeighbors(bb)[ni]->getSecondBead()->getStableIndex();
            beadSet[cumnn + idx] = bindex;
            krep[cumnn+idx] = bb->getRepulsionConst();
            slen[cumnn+idx] = bb->getScreeningLength();
            idx++;

        }
        nneighbors[idb] = idx;
        bubbleSet[idb] = bb->getBead()->getStableIndex();
        radius[idb] = bb->getRadius();
        cumnn+=idx;
//        nintvec[idb] = cumnn;
        idb++;
    }

//    delete [] nintvec;


}

template <class BRepulsionInteractionType>
void BubbleCylinderRepulsion<BRepulsionInteractionType>::deallocate() {
    delete [] beadSet;
    delete [] bubbleSet;
    delete [] krep;
    delete [] slen;
    delete [] nneighbors;
}

template <class BRepulsionInteractionType>
floatingpoint BubbleCylinderRepulsion<BRepulsionInteractionType>::computeEnergy(floatingpoint* coord, bool stretched) {
    
    return _FFType.energy(coord, beadSet, bubbleSet, krep, slen, radius, nneighbors);

}

template <class BRepulsionInteractionType>
void BubbleCylinderRepulsion<BRepulsionInteractionType>::computeForces(floatingpoint *coord, floatingpoint *f) {

_FFType.forces(coord, f, beadSet, bubbleSet, krep, slen, radius, nneighbors);

}

template <class BRepulsionInteractionType>
void BubbleCylinderRepulsion<BRepulsionInteractionType>::computeLoadForces() {
    
    for (auto bb : Bubble::getBubbles()) {
        
        //total number of neighbor cylinders
        int cmax = _neighborList->getNeighbors(bb).size();
        for(int ni = 0; ni < cmax; ni++){
            
            //if part of an MTOC, skip
            if(bb->isMTOC()) {
                
                auto mtoc = (MTOC*)bb->getParent();
                auto filaments = mtoc->getFilaments();
                
                auto f = (Filament*) _neighborList->getNeighbors(bb)[ni]->getParent();
                
                if(find(filaments.begin(), filaments.end(), f) != filaments.end())
                continue;
            }
            
            floatingpoint kRep = bb->getRepulsionConst();
            floatingpoint screenLength = bb->getScreeningLength();
            
            floatingpoint radius = bb->getRadius();
            
            Bead* bd1 = bb->getBead();
            
            //potential acts on second bead unless this is a minus end
            Bead* bd2;
            Bead* bo;
            if(_neighborList->getNeighbors(bb)[ni]->isPlusEnd()) {
                bd2 = _neighborList->getNeighbors(bb)[ni]->getSecondBead();
                bo = _neighborList->getNeighbors(bb)[ni]->getFirstBead();
                
                ///this normal is in the direction of polymerization
                auto normal = normalizeVector(twoPointDirection(bo->vcoordinate(), bd2->vcoordinate()));
                
                //array of coordinate values to update
                auto monSize = SysParams::Geometry().monomerSize[bd2->getType()];
                auto cylSize = SysParams::Geometry().cylinderNumMon[bd2->getType()];
                
                bd2->lfip = 0;
                for (int i = 0; i < cylSize; i++) {
                    
                    auto newCoord = vector<floatingpoint>{bd2->vcoordinate()[0] + i * normal[0] * monSize,
                        bd2->vcoordinate()[1] + i * normal[1] * monSize,
                        bd2->vcoordinate()[2] + i * normal[2] * monSize};
                    
                    // Projection magnitude ratio on the direction of the cylinder
                    // (Effective monomer size) = (monomer size) * proj
                    floatingpoint proj = dotProduct(twoPointDirection(newCoord, bd1->vcoordinate()), normal);
                    if(proj < 0.0) proj = 0.0;
                    
                    floatingpoint loadForce = _FFType.loadForces(bd1, bd2, radius, kRep, screenLength);
                    bd2->loadForcesP[bd2->lfip++] += proj * loadForce;
                }
                //reset lfi
                bd2->lfip = 0;
                
            }
            if(_neighborList->getNeighbors(bb)[ni]->isMinusEnd()) {
                bd2 = _neighborList->getNeighbors(bb)[ni]->getFirstBead();
                bo = _neighborList->getNeighbors(bb)[ni]->getSecondBead();
                
                ///this normal is in the direction of polymerization
                auto normal = normalizeVector(twoPointDirection(bo->vcoordinate(), bd2->vcoordinate()));
                
                //array of coordinate values to update
                auto monSize = SysParams::Geometry().monomerSize[bd2->getType()];
                auto cylSize = SysParams::Geometry().cylinderNumMon[bd2->getType()];
                
                bd2->lfim = 0;
                for (int i = 0; i < cylSize; i++) {
                    
                    auto newCoord = vector<floatingpoint>{bd2->vcoordinate()[0] + i * normal[0] * monSize,
                        bd2->vcoordinate()[1] + i * normal[1] * monSize,
                        bd2->vcoordinate()[2] + i * normal[2] * monSize};
                    
                    // Projection magnitude ratio on the direction of the cylinder
                    // (Effective monomer size) = (monomer size) * proj
                    floatingpoint proj = dotProduct(twoPointDirection(newCoord, bd1->vcoordinate()), normal);
                    if(proj < 0.0) proj = 0.0;
                    
                    floatingpoint loadForce = _FFType.loadForces(bd1, bd2, radius, kRep, screenLength);
                    bd2->loadForcesM[bd2->lfim++] += proj * loadForce;
                }
                //reset lfi
                bd2->lfim = 0;
            }
        }
    }
}


///Template specializations
template floatingpoint BubbleCylinderRepulsion<BubbleCylinderRepulsionExp>::computeEnergy(floatingpoint *coord, bool stretched);
template void BubbleCylinderRepulsion<BubbleCylinderRepulsionExp>::computeForces(floatingpoint *coord, floatingpoint *f);
//template void BubbleCylinderRepulsion<BubbleCylinderRepulsionExp>::computeForcesAux(double *coord, double *f);
template void BubbleCylinderRepulsion<BubbleCylinderRepulsionExp>::computeLoadForces();
template void BubbleCylinderRepulsion<BubbleCylinderRepulsionExp>::vectorize();
template void BubbleCylinderRepulsion<BubbleCylinderRepulsionExp>::deallocate();



