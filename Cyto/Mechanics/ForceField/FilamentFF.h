
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

#ifndef M3SYM_FilamentFF_h
#define M3SYM_FilamentFF_h

#include <vector>

#include "common.h"

#include "ForceField.h"

//FORWARD DECLARATIONS
class FilamentInteractions;

/// Filament FF is an implementation of the [ForceField](@ref ForceField) class that calculates [Filament] (@ref Filament)
/// stretching, bending, and twisting.
class FilamentFF : public ForceField {
 
private:
    vector<unique_ptr<FilamentInteractions>> _filamentInteractionVector; ///< Vector of initialized filament interactions
    
public:
    /// Constructor, intializes stretching, bending, and twisting forces
    FilamentFF(string& stretching, string& bending, string& twisting);
    
    virtual double computeEnergy(double d);
    virtual void computeForces();
    virtual void computeForcesAux();
};

#endif
