
//------------------------------------------------------------------
//  **M3SYM** - Simulation Package for the Mechanochemical
//              Dynamics of Active Networks, 3rd Generation
//
//  Copyright (2015)  Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the MEDYAN web page for more information:
//  http://www.medyan.org
//------------------------------------------------------------------

#ifndef M3SYM_BubbleInteractions_h
#define M3SYM_BubbleInteractions_h

#include "common.h"

//FORWARD DECLARATIONS
class NeighborList;
class Component;
class Bubble;

/// Represents a Bubble interaction with a Bead
class BubbleInteractions {
    
friend class BubbleFF;
    
protected:
    //@{
    /// In the case of an error
    Bubble* _bubbleCulprit;
    Component* _otherCulprit;
    //@}
    
public:
    /// Compute energy of this interaction
    virtual double computeEnergy(double d) = 0;
    /// Compute forces of this interaction
    virtual void computeForces() = 0;
    /// Compute auxiliary forces of this interaction
    virtual void computeForcesAux() = 0;
    
    /// Get the neighbor list for this interaction
    virtual NeighborList* getNeighborList() = 0;
    
    /// Get the name of this interaction
    virtual const string getName() = 0;
};

#endif
