
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v3.0
//
//  Copyright (2015)  Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the MEDYAN web page for more information:
//  http://www.medyan.org
//------------------------------------------------------------------

#ifndef MEDYAN_LinkerInteractions_h
#define MEDYAN_LinkerInteractions_h

#include "common.h"

//FORWARD DECLARATIONS
class Linker;

/// Represents an internal Linker interaction
class LinkerInteractions {
    
friend class LinkerFF;
    
protected:
    /// The linker in the case of an error
    Linker* _linkerCulprit = nullptr;
    
public:
    /// Compute the energy of this interaction
    virtual double computeEnergy(double d) = 0;
    /// Compute the forces of this interaction
    virtual void computeForces() = 0;
    /// Compute the auxiliary forces of this interaction
    virtual void computeForcesAux() = 0;
    
    /// Get the name of this interaction
    virtual const string getName() = 0;
    
};

#endif
