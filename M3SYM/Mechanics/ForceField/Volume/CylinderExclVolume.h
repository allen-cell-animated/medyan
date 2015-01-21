
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

#ifndef M3SYM_CylinderExclVolume_h
#define M3SYM_CylinderExclVolume_h

#include "common.h"

#include "CylinderVolumeInteractions.h"

//FORWARD DECLARATIONS
class Cylinder;

/// Represents an excuded volume interaction between two [Cylinders](@ref Cylinder).
template <class CVolumeInteractionType>
class CylinderExclVolume : public CylinderVolumeInteractions {
    
private:
    CVolumeInteractionType _FFType;
    
public:
    virtual double computeEnergy(Cylinder*, Cylinder*, double d);
    virtual void computeForces(Cylinder*, Cylinder*);
    virtual void computeForcesAux(Cylinder*, Cylinder*);
};

#endif
