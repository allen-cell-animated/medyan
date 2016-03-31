
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

#ifndef M3SYM_LinkerInteractions_h
#define M3SYM_LinkerInteractions_h

#include "common.h"

//FORWARD DECLARATIONS
class Linker;

/// Represents an internal Linker interaction
class LinkerInteractions {
    
friend class LinkerFF;
    
protected:
    /// The linker in the case of an error
    Linker* _linkerCulprit;
    
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
