
//------------------------------------------------------------------
//  **M3SYM** - Simulation Package for the Mechanochemical
//              Dynamics of Active Networks, 3rd Generation
//
//  Copyright (2014) Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the Papoian lab page for more information:
//  http://papoian.chem.umd.edu/
//------------------------------------------------------------------

#ifndef M3SYM_CylinderExclVolume_h
#define M3SYM_CylinderExclVolume_h

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
