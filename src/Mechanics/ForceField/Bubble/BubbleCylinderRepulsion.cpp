
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

#include "BubbleCylinderRepulsion.h"

#include <algorithm> // max

#include "BubbleCylinderRepulsionExp.h"

#include "MTOC.h"
#include "AFM.h"
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
            //if part of an AFM, skip
            else if(bb->isAFM()){
                auto afm = (AFM*)bb->getParent();
                auto filaments = afm->getFilaments();
                
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
            //if part of an AFM, skip
            else if(bb->isAFM()) {
                
                auto afm = (AFM*)bb->getParent();
                auto filaments = afm->getFilaments();
                
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

namespace {

template< typename InteractionType >
void bubbleCylinderRepulsionLoadForce(
    const InteractionType& interaction, floatingpoint radius, floatingpoint kRep, floatingpoint screenLen,
    const Bead& bo, Bead& bd, typename BubbleCylinderRepulsion< InteractionType >::LoadForceEnd end,
    Bead* bbb
) {
    using LoadForceEnd = typename BubbleCylinderRepulsion< InteractionType >::LoadForceEnd;

    auto& loadForces = (end == LoadForceEnd::Plus ? bd.loadForcesP : bd.loadForcesM);
    auto& lfi        = (end == LoadForceEnd::Plus ? bd.lfip        : bd.lfim       );

    // Direction of polymerization
    const auto dir = normalizedVector(bd.coordinate() - bo.coordinate());

    // Array of coordinate values to update
    const auto monSize = SysParams::Geometry().monomerSize   [bd.getType()];
    const auto cylSize = SysParams::Geometry().cylinderNumMon[bd.getType()];

    for (int i = 0; i < cylSize; i++) {

        const auto newCoord = bd.coordinate() + (i * monSize) * dir;

        // Projection magnitude ratio on the direction of the cylinder
        // (Effective monomer size) = (monomer size) * proj
        const auto proj = std::max< floatingpoint >(dot(normalizedVector(bbb->coordinate() - newCoord), dir), 0.0);
        const auto loadForce = interaction.loadForces(bbb, &bd, radius, kRep, screenLen);

        // The load force stored in bead also considers effective monomer size.
        loadForces[i] += proj * loadForce;
    }

    //reset lfi
    lfi = 0;

} // void bubbleCylinderRepulsionLoadForce(...)

} // namespace (anonymous)

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
            //if part of an AFM, skip
            else if(bb->isAFM()) {
                
                auto afm = (AFM*)bb->getParent();
                auto filaments = afm->getFilaments();
                
                auto f = (Filament*) _neighborList->getNeighbors(bb)[ni]->getParent();
                
                if(find(filaments.begin(), filaments.end(), f) != filaments.end())
                    continue;
            }
            
            floatingpoint kRep = bb->getRepulsionConst();
            floatingpoint screenLength = bb->getScreeningLength();
            
            floatingpoint radius = bb->getRadius();
            
            Bead* bd1 = bb->getBead();
            
            Cylinder* c = _neighborList->getNeighbors(bb)[ni];

            if(c->isPlusEnd()) {
                bubbleCylinderRepulsionLoadForce(
                    _FFType, radius, kRep, screenLength,
                    *c->getFirstBead(), *c->getSecondBead(), LoadForceEnd::Plus,
                    bd1
                );
            }
            if(c->isMinusEnd()) {
                bubbleCylinderRepulsionLoadForce(
                    _FFType, radius, kRep, screenLength,
                    *c->getSecondBead(), *c->getFirstBead(), LoadForceEnd::Minus,
                    bd1
                );
            }
        }
    }
}
template< typename InteractionType >
void BubbleCylinderRepulsion< InteractionType >::computeLoadForce(Cylinder* c, LoadForceEnd end) const {
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
            
            Cylinder* cyl = _neighborList->getNeighbors(bb)[ni];
            if(cyl == c) {
                bubbleCylinderRepulsionLoadForce(
                    _FFType, radius, kRep, screenLength,
                    (end == LoadForceEnd::Plus ? *c->getFirstBead() : *c->getSecondBead()),
                    (end == LoadForceEnd::Plus ? *c->getSecondBead() : *c->getFirstBead()),
                    end,
                    bd1
                );
                break;
            }
        } // End loop neighbor of bb
    } // End loop bubbles
}

// Explicit template instantiations
template class BubbleCylinderRepulsion< BubbleCylinderRepulsionExp >;
