
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

#ifndef MEDYAN_CaMKIIingDihedral_h
#define MEDYAN_CaMKIIingDihedral_h

#include "common.h"

#include "CaMKIIingInteractions.h"

//FORWARD DECLARATIONS
class CaMKIIingPoint;

/// Represents an interaction keeping CaMKIIingPoint in dihedral plane
template <class BDihedralInteractionType>
class CaMKIIingDihedral : public CaMKIIingInteractions {
    
private:
    BDihedralInteractionType _FFType;
    
public:
    virtual double computeEnergy(double d);
    virtual void computeForces();
    virtual void computeForcesAux();
    
    virtual const string getName() {return "CaMKIIing Dihedral";}
};

#endif
