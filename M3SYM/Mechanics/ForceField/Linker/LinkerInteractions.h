
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

#ifndef M3SYM_LinkerInteractions_h
#define M3SYM_LinkerInteractions_h

#include "common.h"

//FORWARD DECLARATIONS
class Linker;

/// Represents an internal Linker interaction
class LinkerInteractions {
    
public:
    /// Compute the energy of this interaction
    virtual double computeEnergy(Linker*,  double d) = 0;
    /// Compute the forces of this interaction
    virtual void computeForces(Linker*) = 0;
    /// Compute the auxiliary forces of this interaction
    virtual void computeForcesAux(Linker*) = 0;
    
    /// Get the name of this interaction
    virtual const string getName() = 0;
    
};

#endif
