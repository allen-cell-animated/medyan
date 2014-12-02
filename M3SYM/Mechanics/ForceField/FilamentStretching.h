
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

#ifndef M3SYM_FilamentStretching_h
#define M3SYM_FilamentStretching_h

#include "common.h"

#include "FilamentInteractions.h"

/// Represents a filament stretching interaction
template <class FStretchingInteractionType>
class FilamentStretching : public FilamentInteractions {
    
private:
    FStretchingInteractionType _FFType; 
    
public:
    virtual double computeEnergy(Filament*, double d);
    virtual void computeForces(Filament*);
    virtual void computeForcesAux(Filament*);
};


#endif
