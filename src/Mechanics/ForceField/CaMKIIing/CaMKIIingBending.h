
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

#ifndef MEDYAN_CaMKIIingBending_h
#define MEDYAN_CaMKIIingBending_h

#include "common.h"

#include "CaMKIIingInteractions.h"

//FORWARD DECLARATIONS
class CaMKIIingPoint;

/// Represents an interaction maintaining a CaMKIIingPoint angle (~270 for Arp2/3)
template <class BBendingInteractionType>
class CaMKIIingBending : public CaMKIIingInteractions {
    
private:
    BBendingInteractionType _FFType;
    
public:
    virtual double computeEnergy(double d);
    virtual void computeForces();
    virtual void computeForcesAux();
    
    virtual const string getName() {return "CaMKIIing Bending";}
};

#endif
