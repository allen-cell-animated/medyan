
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

#ifndef MEDYAN_BoundaryCylinderRepulsionExpIn_h
#define MEDYAN_BoundaryCylinderRepulsionExpIn_h

#include <vector>
#include <cmath>

#include "common.h"

//FORWARD DECLARATIONS
class Bead;

/// A exponential repulsive potential used by the BoundaryCylinderRepulsion template.
class BoundaryCylinderRepulsionExpIn {
    
public:
    //TODO needs implementation @{
    floatingpoint energy(floatingpoint *coord, floatingpoint *f, int *beadSet,
                  floatingpoint *krep, floatingpoint *slen, int *nneighbors);

    floatingpoint energy(floatingpoint *coord, floatingpoint *f, int *beadSet,
                  floatingpoint *krep, floatingpoint *slen, int *nnneighbors, floatingpoint d);

    void forces(floatingpoint *coord, floatingpoint *f, int *beadSet,
                floatingpoint *krep, floatingpoint *slen, int *nneighbors);
    //@}
    floatingpoint energy(Bead*, floatingpoint, floatingpoint, floatingpoint);
    void forces(Bead*, floatingpoint, vector<floatingpoint>& norm, floatingpoint, floatingpoint);
    void forcesAux(Bead*, floatingpoint, vector<floatingpoint>& norm, floatingpoint, floatingpoint);
    floatingpoint loadForces(floatingpoint, floatingpoint, floatingpoint);
};

#endif