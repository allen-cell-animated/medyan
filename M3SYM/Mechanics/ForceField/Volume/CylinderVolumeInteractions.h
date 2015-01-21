
//------------------------------------------------------------------
//  **M3SYM** - Simulation Package for the Mechanochemical
//              Dynamics of Active Networks, 3rd Generation
//
//  Copyright (2014) Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the Papoian lab page for installation and documentation:
//  http://papoian.chem.umd.edu/
//------------------------------------------------------------------

#ifndef M3SYM_CylinderVolInteractions_h
#define M3SYM_CylinderVolInteractions_h

#include "common.h"

#include "NeighborListContainer.h"
#include "SystemParameters.h"

//FORWARD DECLARATIONS
class Cylinder;

/// Represents a volume interaction between [Cylinders](@ref Cylinder).
class CylinderVolumeInteractions : public CCNLContainer {
private:
    string _name; ///< Name of interaction
    
public:
    ///Constructor, initializes a cylinder neighbor list
    CylinderVolumeInteractions()
        : CCNLContainer(SystemParameters::Mechanics().VolumeCutoff) {}
    
    /// Compute the energy of this interaction
    virtual double computeEnergy(Cylinder*, Cylinder*,  double d) = 0;
    /// Compute the forces of this interaction
    virtual void computeForces(Cylinder*, Cylinder*) = 0;
    /// Compute the auxiliary forces of this interaction
    virtual void computeForcesAux(Cylinder*, Cylinder*) = 0;
    
    /// Get name of interaction
    const string& getName() {return _name;}
};


#endif
