
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v3.2.1
//
//  Copyright (2015-2018)  Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the MEDYAN web page for more information:
//  http://www.medyan.org
//------------------------------------------------------------------

#ifndef MEDYAN_FilamentInteractions_h
#define MEDYAN_FilamentInteractions_h

#include <iostream>

#include "common.h"

//FORWARD DECLARATIONS
class Filament;

/// Represents an internal Filament interaction
class FilamentInteractions {
    
friend class FilamentFF;
    
public:
    /// The filament in the case of an error
    static Filament* _filamentCulprit;
    
    ///Vectorize the bead interactions for minimization
    virtual void vectorize() = 0;
    ///Deallocate the vectorized data
    virtual void deallocate() = 0;
    
    /// Compute the energy of this interaction
    virtual double computeEnergy(double *coord, double *f, double d) = 0;
    /// Compute the forces of this interaction
    virtual void computeForces(double *coord, double *f) = 0;
    
    /// Get the name of this interaction
    virtual const string getName() = 0;

    /// Get culprit information. Used in CUDA.
//    virtual void whoisculprit();
};

#endif
