
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

#ifndef MEDYAN_FilamentStretching_h
#define MEDYAN_FilamentStretching_h

#include "common.h"

#include "FilamentInteractions.h"

/// Represents a Filament stretching interaction
template <class FStretchingInteractionType>
class FilamentStretching : public FilamentInteractions {
    
private:
    FStretchingInteractionType _FFType; 
    
public:
    virtual double computeEnergy(bool stretched) override;
    virtual void computeForces();
    virtual void computeForcesAux();
    
    virtual const string getName() {return "Filament Stretching";}
};


#endif
