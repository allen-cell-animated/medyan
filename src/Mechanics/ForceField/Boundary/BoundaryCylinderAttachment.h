
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
    
    int *beadSet;
    
    ///Array describing the constants in calculation
    double *kattr;
    double *pins; ///< coordinates of pins for each bead
    
public:
    
    ///Array describing indexed set of interactions
    ///For filaments, this is a 1-point potential (only the active bead is updated)
    const static int n = 1;
    
    virtual void vectorize();
    virtual void deallocate();
    
    virtual double computeEnergy(double *coord) override;
    //@{
    /// Tepulsive force calculation
    virtual void computeForces(double *coord, double *f);
    virtual void computeLoadForces() {};
    //@}
    
    /// Get the neighbor list for this interaction
    virtual NeighborList* getNeighborList() {return nullptr;}
    
    virtual const string getName() {return "Boundary-Cylinder Attachment";}
};
#endif
