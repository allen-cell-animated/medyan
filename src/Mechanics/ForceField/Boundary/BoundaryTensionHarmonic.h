/*
 * BoundaryTensionHarmonic.h
 *
 *  Created on: Jul 19, 2016
 *      Author: james
 */

#ifndef MECHANICS_FORCEFIELD_BOUNDARY_BOUNDARYTENSIONHARMONIC_H_
#define MECHANICS_FORCEFIELD_BOUNDARY_BOUNDARYTENSIONHARMONIC_H_

#include "common.h"

#include "BoundaryElementImpl.h"

/// A harmonic tension used by Boundary
class BoundaryTensionHarmonic {

public:
    double energy(double, double);
    double energy(double, double, double);

    void forces(SphereBoundaryElement*, double, double);
    void forcesAux(SphereBoundaryElement*, double, double);
};



#endif /* MECHANICS_FORCEFIELD_BOUNDARY_BOUNDARYTENSIONHARMONIC_H_ */
