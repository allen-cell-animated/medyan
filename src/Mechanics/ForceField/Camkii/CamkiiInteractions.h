
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

#ifndef MEDYAN_CamkiiInteractions_h
#define MEDYAN_CamkiiInteractions_h

#include <iostream>

#include "common.h"

//FORWARD DECLARATIONS
class Filament;
class Camkii;

/// Represents an internal Filament interaction
class CamkiiInteractions {
    
friend class CamkiiFF;
    
protected:
    /// The filament in the case of an error
    Camkii* _camkiiCulprit = nullptr;

public:
    /// Compute the energy of this interaction
    virtual double computeEnergy(double d) = 0;
    /// Compute forces of this interaction
    virtual void computeForces() = 0;
    /// Compute auxiliary forces of this interaction
    virtual void computeForcesAux() = 0;
    
    /// Get the name of this interaction
    virtual const string getName() = 0;
};

#endif
