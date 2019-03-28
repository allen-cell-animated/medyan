
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

#ifndef MEDYAN_BubbleInteractions_h
#define MEDYAN_BubbleInteractions_h

#include "common.h"

//FORWARD DECLARATIONS
class NeighborList;
class Component;
class Bubble;

/// Represents a Bubble interaction with a Bead
class BubbleInteractions {
    
friend class BubbleFF;
    
//protected:
//    //@{
//    /// In the case of an error
//    Bubble* _bubbleCulprit = nullptr;
//    Component* _otherCulprit = nullptr;
//    //@}
    
public:
    
    //@{
    /// In the case of an error
    static Bubble* _bubbleCulprit;
    static Component* _otherCulprit;
    //@}
    
    ///Vectorize the bead interactions for minimization
    virtual void vectorize() = 0;
    ///Deallocate the vectorized data
    virtual void deallocate() = 0;
    
    /// Compute energy of this interaction
    virtual double computeEnergy(double *coord, bool stretched) = 0;
    /// Compute forces of this interaction
    virtual void computeForces(double *coord, double *f) = 0;
    /// Compute auxiliary forces of this interaction
    //virtual void computeForcesAux() = 0;
    
    /// Compute the load forces on beads from this interaction
    virtual void computeLoadForces() = 0;
    
    /// Get the neighbor list for this interaction
    virtual NeighborList* getNeighborList() = 0;
    
    /// Get the name of this interaction
    virtual const string getName() = 0;
};

#endif
