
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

#ifndef MEDYAN_TriangleBeadExclVolume_h
#define MEDYAN_TriangleBeadExclVolume_h

#include "common.h"

#include "TriangleBeadVolumeInteractions.h"
#include "NeighborListImpl.h"

#include "SysParams.h"

//FORWARD DECLARATIONS
class Triangle;
class Bead;

/// Represents an excuded volume interaction between two [Cylinders](@ref Cylinder).
template <class TriangleBeadExclVolumeInteractionType>
class TriangleBeadExclVolume : public TriangleBeadVolumeInteractions {
    
private:
    TriangleBeadExclVolumeInteractionType _FFType;
    CylinderCylinderNL* _neighborList;  ///< TODO: Change neighbor list.
    
public:
    ///Constructor
    TriangleBeadExclVolume() {
        //v TODO: change this.
        _neighborList = new CylinderCylinderNL(SysParams::Mechanics().VolumeCutoff);
    }
    
    virtual double computeEnergy(double d);
    virtual void computeForces();
    virtual void computeForcesAux();

    /// Get the neighbor list for this interaction
    virtual NeighborList* getNeighborList() { return _neighborList; }
    
    virtual const string getName() {return "Triangle Bead Excluded Volume";}
};

#endif
