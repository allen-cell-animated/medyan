
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

#ifndef MEDYAN_CylinderExclVolume_h
#define MEDYAN_CylinderExclVolume_h

#include "common.h"

#include "CylinderVolumeInteractions.h"
#include "NeighborListImpl.h"

#include "SysParams.h"

//FORWARD DECLARATIONS
class Cylinder;

/// Represents an excuded volume interaction between two [Cylinders](@ref Cylinder).
template <class CVolumeInteractionType>
class CylinderExclVolume : public CylinderVolumeInteractions {
    
private:
    CVolumeInteractionType _FFType;
    CylinderCylinderNL* _neighborList;  ///< Neighbor list of cylinders
    
public:
    ///Constructor
    CylinderExclVolume() {
        _neighborList = new CylinderCylinderNL(SysParams::Mechanics().VolumeCutoff);
    }
    
    virtual double computeEnergy(double d);
    virtual void computeForces();
    virtual void computeForcesAux();

    /// Get the neighbor list for this interaction
    virtual NeighborList* getNeighborList() {return _neighborList;}
    
    virtual const string getName() {return "Cylinder Excluded Volume";}
};

#endif
