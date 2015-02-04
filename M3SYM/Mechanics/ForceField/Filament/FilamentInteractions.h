
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

#ifndef M3SYM_FilamentInteractions_h
#define M3SYM_FilamentInteractions_h

#include <iostream>

#include "common.h"

//FORWARD DECLARATIONS
class Filament;

/// Represents an internal Filament interaction
class FilamentInteractions {
private:
    string _name;

public:
    /// Compute the energy of this interaction
    virtual double computeEnergy(Filament*,  double d) = 0;
    /// Compute forces of this interaction
    virtual void computeForces(Filament*) = 0;
    /// Compute auxiliary forces of this interaction
    virtual void computeForcesAux(Filament*) = 0;
    
    /// Get the name of this interaction
    const string& getName() {return _name;}
    
};

#endif
