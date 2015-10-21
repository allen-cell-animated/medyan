
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

#ifndef M3SYM_BubbleBubbleRepulsion_h
#define M3SYM_BubbleBubbleRepulsion_h

#include <vector>

#include "common.h"

#include "BubbleInteractions.h"
#include "NeighborListImpl.h"

#include "SysParams.h"

//FORWARD DECLARATIONS
class Bead;

/// Represents a repulsive interaction between two Bubbles.
template <class BRepulsionInteractionType>
class BubbleBubbleRepulsion : public BubbleInteractions {
    
private:
    BRepulsionInteractionType _FFType;
    BubbleBubbleNL* _neighborList; ///<Neighbor list of Bubble-Bubble
public:
    
    /// Constructor
    BubbleBubbleRepulsion() {
        _neighborList = new BubbleBubbleNL(SysParams::Mechanics().BubbleCutoff);
    }
    
    virtual double computeEnergy(double d);
    
    virtual void computeForces();
    virtual void computeForcesAux();
    
    /// Get the neighbor list for this interaction
    virtual NeighborList* getNeighborList() {return _neighborList;}
    
    virtual const string getName() {return "Bubble-Bubble Repulsion";}
};
#endif
