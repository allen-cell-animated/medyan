
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

#ifndef MEDYAN_AFMAttachment_h
#define MEDYAN_AFMAttachment_h

#include <vector>

#include "common.h"

#include "BubbleInteractions.h"

#include "SysParams.h"

//FORWARD DECLARATIONS
class AFM;
class Bead;

/// Represents an attachment potential of a AFM.
template <class AFMInteractionType>
class AFMAttachment : public BubbleInteractions {
    
private:
    AFMInteractionType _FFType;

    int numInteractions_ = 0;
    std::vector< int > beadSet_;
    ///Array describing the constants in calculation
    std::vector< floatingpoint > radii_;
    std::vector< floatingpoint > kstr_;
    floatingpoint *pos1;
    floatingpoint *pos2;
    
    
public:
    
    ///Array describing indexed set of interactions
    ///For AFM, this is a 1-bead potential
    const static int n = 1;
    
    virtual void vectorize(const FFCoordinateStartingIndex&) override;
    virtual void deallocate() {}
    
    virtual floatingpoint computeEnergy(floatingpoint *coord, bool stretched) override;
    virtual void computeForces(floatingpoint *coord, floatingpoint *f);
    //virtual void computeForcesAux(double *coord, double *f);
    
    virtual void computeLoadForces() {return;}
    
    /// Get the neighbor list for this interaction
    virtual NeighborList* getNeighborList() {return nullptr;}
    
    virtual const string getName() {return "AFM Attachment";}
};

#endif


