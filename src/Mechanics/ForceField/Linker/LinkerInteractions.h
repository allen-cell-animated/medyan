
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

#ifndef MEDYAN_LinkerInteractions_h
#define MEDYAN_LinkerInteractions_h

#include "common.h"

//FORWARD DECLARATIONS
class Linker;

/// Represents an internal Linker interaction
class LinkerInteractions {
    
friend class LinkerFF;
    
public:
    /// The linker in the case of an error
    static Linker* _linkerCulprit;
    
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

    virtual void assignforcemags(){};
    
};

#endif
