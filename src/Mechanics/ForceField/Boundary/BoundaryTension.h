/*
 * BoundaryTension.h
 *
 *  Created on: Jul 19, 2016
 *      Author: jl135
 */

#ifndef MECHANICS_FORCEFIELD_BOUNDARY_BOUNDARYTENSION_H_
#define MECHANICS_FORCEFIELD_BOUNDARY_BOUNDARYTENSION_H_

#include "common.h"

#include "BoundaryTensionHarmonic.h"

#include "BoundaryInteractions.h"

/// Template for type of boundary tension (harmonic, etc.)
template <class BTensionType>
class BoundaryTension : public BoundaryInteractions {

private:
	BTensionType _FFType;

public:
    virtual double computeEnergy(double d);
    virtual void computeForces();
    virtual void computeForcesAux();

    virtual const string getName() {return "Tension Boundary is being Applied";}
};


#endif /* MECHANICS_FORCEFIELD_BOUNDARY_BOUNDARYTENSION_H_ */
