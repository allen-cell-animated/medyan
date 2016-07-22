
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v3.0
//
//  Copyright (2015)  Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the MEDYAN web page for more information:
//  http://www.medyan.org
//------------------------------------------------------------------

#ifndef MEDYAN_BoundaryCylinderAttachment_h
#define MEDYAN_BoundaryCylinderAttachment_h

#include <vector>

#include "common.h"

#include "BoundaryInteractions.h"

#include "SysParams.h"

//FORWARD DECLARATIONS
class BoundaryElement;
class Bead;

/// Represents an attractive interaction between a cylinder and its pin point near a boundary.
/// @note - This BoundaryInteractions implementation really does not involve Boundaries at all.
///         The pinned position of the beads is preset by the special protocol in which Bead
///         elements are chosen with the criteria of being a certain distance away from the Boundary.
///         In reality, any Bead could be pinned at a certain position, but for now associating
///         this potential with a Boundary makes more sense.

template <class BAttachmentInteractionType>
class BoundaryCylinderAttachment : public BoundaryInteractions {
    
private:
    BAttachmentInteractionType _FFType;
public:
    
    /// Constructor
    BoundaryCylinderAttachment() {}
    
    virtual double computeEnergy(double d);
    
    virtual void computeForces();
    virtual void computeForcesAux();
    
    virtual void computeLoadForces() {}
    
    /// Just return null - no neighbor list here
    virtual NeighborList* getNeighborList() {return nullptr;}
    
    virtual const string getName() {return "Boundary-Cylinder Attachment";}
};
#endif
