
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

#ifndef MEDYAN_TriangleCylinderExclVolume_h
#define MEDYAN_TriangleCylinderExclVolume_h

#include "common.h"

#include "TriangleCylinderVolumeInteractions.h"
#include "NeighborListImpl.h"

#include "SysParams.h"

//FORWARD DECLARATIONS
class Triangle;
class Cylinder;
class Bead;

/// Represents an excuded volume interaction between a triangle and a cylinder (bead).
template <class TriangleCylinderExclVolumeInteractionType>
class TriangleCylinderExclVolume : public TriangleCylinderVolumeInteractions {
    
private:
    TriangleCylinderExclVolumeInteractionType _FFType;
    TriangleCylinderNL* _neighborList;  ///< Neighbor list of cylinders
    
public:
    ///Constructor
    TriangleCylinderExclVolume() {
        _neighborList = new TriangleCylinderNL(SysParams::Mechanics().MemCylinderVolumeCutoff);
    }
    
    virtual double computeEnergy(double d);
    virtual void computeForces();
    virtual void computeForcesAux();

    /// Get the neighbor list for this interaction
    virtual NeighborList* getNeighborList() { return _neighborList; }
    
    virtual const string getName() {return "Triangle Cylinder Excluded Volume";}
};

#endif
