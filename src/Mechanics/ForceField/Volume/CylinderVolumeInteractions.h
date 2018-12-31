
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

#ifndef MEDYAN_CylinderVolumeInteractions_h
#define MEDYAN_CylinderVolumeInteractions_h

#include "common.h"

//FORWARD DECLARATIONS
class NeighborList;
class Cylinder;

/// Represents a volume interaction between [Cylinders](@ref Cylinder).
class CylinderVolumeInteractions {

friend class CylinderVolumeFF;
    
protected:
    //@{
    /// The cylinder culprits in the case of an error
    Cylinder* _cylinderCulprit1 = nullptr;
    Cylinder* _cylinderCulprit2 = nullptr;
    //@}
    
public:
    /// Compute the energy of this interaction
    virtual double computeEnergy(bool stretched) = 0;
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
