
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

#ifndef MEDYAN_MTOCAttachment_h
#define MEDYAN_MTOCAttachment_h

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
    virtual floatingpoint computeEnergy(floatingpoint d);
    
    virtual void computeForces();
    virtual void computeForcesAux();
    
    virtual void computeLoadForces() {return;}
    
    /// Get the neighbor list for this interaction
    virtual NeighborList* getNeighborList() {return nullptr;}
    
    virtual const string getName() {return "MTOC Attachment";}
};

#endif
