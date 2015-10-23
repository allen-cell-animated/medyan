
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

#ifndef M3SYM_MTOCAttachment_h
#define M3SYM_MTOCAttachment_h

#include <vector>

#include "common.h"

#include "BubbleInteractions.h"

#include "SysParams.h"

//FORWARD DECLARATIONS
class MTOC;
class Bead;

/// Represents an attachment potential of a MTOC.
template <class MTOCInteractionType>
class MTOCAttachment : public BubbleInteractions {
    
private:
    MTOCInteractionType _FFType;
public:
    virtual double computeEnergy(double d);
    
    virtual void computeForces();
    virtual void computeForcesAux();
    
    /// Get the neighbor list for this interaction
    virtual NeighborList* getNeighborList() {return nullptr;}
    
    virtual const string getName() {return "MTOC Attachment";}
};

#endif
