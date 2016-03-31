
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

#ifndef M3SYM_FilamentStretching_h
#define M3SYM_FilamentStretching_h

#include "common.h"

#include "FilamentInteractions.h"

/// Represents a Filament stretching interaction
template <class FStretchingInteractionType>
class FilamentStretching : public FilamentInteractions {
    
private:
    FStretchingInteractionType _FFType; 
    
public:
    virtual double computeEnergy(double d);
    virtual void computeForces();
    virtual void computeForcesAux();
    
    virtual const string getName() {return "Filament Stretching";}
};


#endif
