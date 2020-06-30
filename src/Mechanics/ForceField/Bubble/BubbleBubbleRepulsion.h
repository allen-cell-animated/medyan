
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

#ifndef MEDYAN_BubbleBubbleRepulsion_h
#define MEDYAN_BubbleBubbleRepulsion_h

#include <vector>

#include "common.h"

#include "BubbleInteractions.h"
#include "NeighborListImpl.h"

#include "SysParams.h"

//FORWARD DECLARATIONS
class Bead;

/// Represents a repulsive interaction between two [Bubbles](@ref Bubble).
template <class BRepulsionInteractionType>
class BubbleBubbleRepulsion : public BubbleInteractions {
    
private:
    BRepulsionInteractionType _FFType;
    BubbleBubbleNL* _neighborList; ///<Neighbor list of Bubble-Bubble

    std::size_t beadStartIdx_; // The starting index of bead in the vectorization

public:
    
    /// Constructor
    BubbleBubbleRepulsion() {
        _neighborList = new BubbleBubbleNL(SysParams::Mechanics().BubbleCutoff);
    }
    
    virtual void vectorize(const FFCoordinateStartingIndex&) override;
    virtual void deallocate();
    
    virtual floatingpoint computeEnergy(floatingpoint *coord, bool stretched) override;
    virtual void computeForces(floatingpoint *coord, floatingpoint *f);
    //virtual void computeForcesAux(double *coord, double *f);
    
    virtual void computeLoadForces() {return;}
    
    /// Get the neighbor list for this interaction
    virtual NeighborList* getNeighborList() {return _neighborList;}
    
    virtual const string getName() {return "Bubble-Bubble Repulsion";}
};
#endif
