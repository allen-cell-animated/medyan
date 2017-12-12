
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

#ifndef MEDYAN_TriangleBeadVolumeInteractions_h
#define MEDYAN_TriangleBeadVolumeInteractions_h

#include "common.h"

//FORWARD DECLARATIONS
class NeighborList;
class Triangle;
class Bead;

/// Represents a volume interaction between [Cylinders](@ref Cylinder).
class TriangleBeadVolumeInteractions {

friend class TriangleBeadVolumeFF;
    
protected:
    //@{
    /// The triangle and bead culprits in the case of an error
    Triangle* _triangleCulprit = nullptr;
    Bead* _beadCulprit = nullptr;
    //@}
    
public:
    /// Compute the energy of this interaction
    virtual double computeEnergy(double d) = 0;
    /// Compute the forces of this interaction
    virtual void computeForces() = 0;
    /// Compute the auxiliary forces of this interaction
    virtual void computeForcesAux() = 0;
    
    /// Get the neighbor list for this interaction
    virtual NeighborList* getNeighborList() = 0;
    
    /// Get the name of this interaction
    virtual const string getName() = 0;
};


#endif
