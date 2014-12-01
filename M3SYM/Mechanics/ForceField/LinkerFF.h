
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

#ifndef M3SYM_LinkerFF_h
#define M3SYM_LinkerFF_h

#include <vector>

#include "common.h"

#include "ForceField.h"

//FORWARD DECLARATIONS
class LinkerInteractions;

/// Linker FF is an implementation of the [ForceField](@ref ForceField) class that calculates [Linker] (@ref Linker)
/// stretching, bending, and twisting.
class LinkerFF : public ForceField {
    
private:
    vector<unique_ptr<LinkerInteractions>> _linkerInteractionVector; ///< Vector of initialized linker interactions
    
public:
    /// Constructor, intializes stretching, bending, and twisting forces
    LinkerFF(string& stretching, string& bending, string& twisting );
    
    virtual double computeEnergy(double d);
    virtual void computeForces();
    virtual void computeForcesAux();
    
};

#endif