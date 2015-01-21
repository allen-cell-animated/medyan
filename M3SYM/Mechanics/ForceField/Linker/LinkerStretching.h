
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

#ifndef M3SYM_LinkerStretching_h
#define M3SYM_LinkerStretching_h

#include "common.h"

#include "LinkerInteractions.h"

//FORWARD DECLARATIONS
class Linker;

/// Represents a Linker stretching interaction
template <class LStretchingInteractionType>
class LinkerStretching : public LinkerInteractions {
    
private:
    LStretchingInteractionType _FFType;
    
public:
    virtual double computeEnergy(Linker*, double d);
    virtual void computeForces(Linker*);
    virtual void computeForcesAux(Linker*);
};

#endif
