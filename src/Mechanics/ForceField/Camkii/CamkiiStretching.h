#ifdef CAMKII
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v3.1
//
//  Copyright (2015-2016)  Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the MEDYAN web page for more information:
//  http://www.medyan.org
//------------------------------------------------------------------

#ifndef MEDYAN_CamkiiStretching_h
#define MEDYAN_CamkiiStretching_h

#include "common.h"

#include "CamkiiInteractions.h"

/// Represents a Camkii stretching interaction
template <class FStretchingInteractionType>
class CamkiiStretching : public CamkiiInteractions {


private:
    //TODO
    FStretchingInteractionType _FFType; 
    
public:
    virtual double computeEnergy(double d);
    virtual void computeForces();
    virtual void computeForcesAux();
    
    virtual const string getName() {return "Camkii Stretching";}
};

#endif
#endif
