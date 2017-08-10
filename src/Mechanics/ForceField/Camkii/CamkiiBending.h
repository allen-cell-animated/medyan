
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

#ifndef MEDYAN_CamkiiBending_h
#define MEDYAN_CamkiiBending_h

#include "common.h"

#include "CamkiiInteractions.h"

//FORWARD DECLARATIONS
class Camkii;

/// Represents a Camkii bending interaction
template <class FBendingInteractionType>
class CamkiiBending : public CamkiiInteractions {
    
private:
    FBendingInteractionType _FFType;
    
public:
    virtual double computeEnergy(double d);
    virtual void computeForces();
    virtual void computeForcesAux();
    
    virtual const string getName() {return "Camkii Bending";}
};

#endif
