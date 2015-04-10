
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

#ifndef M3SYM_FilamentBending_h
#define M3SYM_FilamentBending_h

#include "common.h"

#include "FilamentInteractions.h"

//FORWARD DECLARATIONS
class Filament;

/// Represents a Filament bending interaction
template <class FBendingInteractionType>
class FilamentBending : public FilamentInteractions {
    
private:
    FBendingInteractionType _FFType;
    
public:
    virtual double computeEnergy(Filament*, double d);
    virtual void computeForces( Filament*);
    virtual void computeForcesAux( Filament*);
    
    virtual const string getName() {return "Bending";}
};

#endif
