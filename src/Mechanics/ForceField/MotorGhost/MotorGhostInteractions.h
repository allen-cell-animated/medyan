
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v4.0
//
//  Copyright (2015-2018)  Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the MEDYAN web page for more information:
//  http://www.medyan.org
//------------------------------------------------------------------


#ifndef MEDYAN_MotorGhostInteractions_h
#define MEDYAN_MotorGhostInteractions_h

#include "common.h"
#include "Mechanics/ForceField/Types.hpp"

//FORWARD DECLARATIONS
class MotorGhost;

/// Represents an internal MotorGhost interaction
class MotorGhostInteractions {
    
friend class MotorGhostFF;
    
public:
    /// The motor ghost in the case of an error
    static MotorGhost* _motorCulprit;

    ///Vectorize the bead interactions for minimization
    virtual void vectorize(const FFCoordinateStartingIndex&) = 0;
    ///Deallocate the vectorized data
    virtual void deallocate() = 0;
    
    /// Compute the energy of this interaction
    virtual floatingpoint computeEnergy(floatingpoint *coord) = 0;
    /// Compute the forces of this interaction
    virtual void computeForces(floatingpoint *coord, floatingpoint *f) = 0;

    /// Get the name of this interaction
    virtual const string getName() = 0;
    /// Assign stretchforces for RateChangerImpl.
    virtual void assignforcemags(){};

    static vector<floatingpoint> individualenergies;
	static vector<floatingpoint> tpdistvec;
	static vector<floatingpoint> eqlvec;
	static vector<floatingpoint> kstrvec;
};

#endif
