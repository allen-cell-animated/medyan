
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
    
    int *beadSet;
    ///Array describing the constants in calculation
    double *kstr;
    double *pos1;
    double *pos2;
    
    
public:
    
    ///Array describing indexed set of interactions
    ///For AFM, this is a 1-bead potential
    const static int n = 1;
    
    virtual void vectorize();
    virtual void deallocate();
    
    virtual double computeEnergy(double *coord, double *f, double d);
    virtual void computeForces(double *coord, double *f);
    //virtual void computeForcesAux(double *coord, double *f);
    
    virtual void computeLoadForces() {return;}
    
    /// Get the neighbor list for this interaction
    virtual NeighborList* getNeighborList() {return nullptr;}
    
    virtual const string getName() {return "AFM Attachment";}
};

#endif


