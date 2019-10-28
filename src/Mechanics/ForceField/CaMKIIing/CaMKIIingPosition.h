
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

#ifndef MEDYAN_CaMKIIingPosition_h
#define MEDYAN_CaMKIIingPosition_h

#include "common.h"

#include "CaMKIIingInteractions.h"

//FORWARD DECLARATIONS
class CaMKIIingPoint;

/// Represents an interaction fixing a Cylinder anchored by a CaMKIIingPoint on the parent.
template <class CStretchingInteractionType>
class CaMKIIingPosition : public CaMKIIingInteractions {
    
private:
    CStretchingInteractionType _FFType;
    
public:
    virtual double computeEnergy(double d);
    virtual void computeForces();
    virtual void computeForcesAux();
    
    virtual const string getName() {return "CaMKIIing Position";}
};

#endif

