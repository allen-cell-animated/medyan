
//------------------------------------------------------------------
//  **M3SYM** - Simulation Package for the Mechanochemical
//              Dynamics of Active Networks, 3rd Generation
//
//  Copyright (2014) Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the Papoian lab page for installation and documentation:
//  http://papoian.chem.umd.edu/
//------------------------------------------------------------------

#ifndef M3SYM_LinkerInteractions_h
#define M3SYM_LinkerInteractions_h

#include "common.h"

//FORWARD DECLARATIONS
class Linker;

/// LinkerInteractions class represents an internal linker interaction
class LinkerInteractions {
private:
    string _name; ///< Name of interaction
    
public:
    /// Compute the energy of this interaction
    virtual double computeEnergy(Linker*,  double d) = 0;
    /// Compute the forces of this interaction
    virtual void computeForces(Linker*) = 0;
    /// Compute the auxiliary forces of this interaction
    virtual void computeForcesAux(Linker*) = 0;
    
    /// Get the name of his interaction
    const string& getName() {return _name;}
    
};

#endif
