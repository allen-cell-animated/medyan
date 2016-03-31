
//------------------------------------------------------------------
//  **M3SYM** - Simulation Package for the Mechanochemical
//              Dynamics of Active Networks, 3rd Generation
//
//  Copyright (2015)  Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the MEDYAN web page for more information:
//  http://www.medyan.org
//------------------------------------------------------------------

#ifndef M3SYM_FilamentInteractions_h
#define M3SYM_FilamentInteractions_h

#include <iostream>

#include "common.h"

//FORWARD DECLARATIONS
class Filament;

/// Represents an internal Filament interaction
class FilamentInteractions {
    
friend class FilamentFF;
    
protected:
    /// The filament in the case of an error
    Filament* _filamentCulprit;

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
