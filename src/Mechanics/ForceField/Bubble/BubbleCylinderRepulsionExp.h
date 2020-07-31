
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v4.0
//
//  Copyright (2015-2018)  Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the MEDYAN web page for more information:
//  http://www.medyan.org
//------------------------------------------------------------------

#ifndef MEDYAN_BubbleCylinderRepulsionExp_h
#define MEDYAN_BubbleCylinderRepulsionExp_h

#include <vector>

#include "common.h"
#include "MathFunctions.h"

//FORWARD DECLARATIONS
class Bead;

/// A exponential repulsive potential used by the BubbleCylinderRepulsion template.
class BubbleCylinderRepulsionExp {
    
public:
    floatingpoint energy(floatingpoint *coord, int *beadSet, int *bubbleSet,
                         floatingpoint *krep, floatingpoint *slen, floatingpoint *radius, int *nneighbors);
    
    [[deprecated]] floatingpoint energy(floatingpoint *coord, floatingpoint *f, int *beadSet, int *bubbleSet,
                         floatingpoint *krep, floatingpoint *slen, floatingpoint *radius, int *nnneighbors, floatingpoint d);
    
    void forces(floatingpoint *coord, floatingpoint *f, int *beadSet, int *bubbleSet,
                floatingpoint *krep, floatingpoint *slen, floatingpoint *radius, int *nneighbors);
//    void forcesAux(Bead*, Bead*, double, double, double);

    template< typename VT1, typename VT2 >
	floatingpoint loadForces(
        const VT1&    coord1,
        const VT2&    coord2,
        floatingpoint radius,
        floatingpoint kRep,
        floatingpoint screenLength
    ) const {

        const auto dist = mathfunc::distance(coord1, coord2);

        floatingpoint effd = dist - radius;

        floatingpoint R = -effd / screenLength;
        return kRep * exp(R) / screenLength;

    }
};

#endif
