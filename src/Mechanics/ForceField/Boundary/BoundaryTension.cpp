/*
 * BoundaryTension.cpp
 *
 *  Created on: Jul 19, 2016
 *      Author: jl135
 */
#include "BoundaryTension.h"

#include "BoundaryTensionHarmonic.h"

#include "BoundaryInteractions.h"

#include "BoundaryElement.h"


template <class BTensionType>
double BoundaryTension<BTensionType>::computeEnergy(double d) {

    double U = 0;
    double U_i;

    for (auto be: BoundaryElement::getBoundaryElements()){
    //for(BoundaryElement* be: BoundaryElement::getBoundaryElements()){ }
    // vector of boundary elements = BoundaryElement::getBoundaryElements()
    	    //	for(BoundaryElement* be : vector) {
    	    //	}

    	SphereBoundaryElement* sb = (SphereBoundaryElement*)be;
    	U_i = 0;
    	if (d == 0.0){

        		double d_radius = sb->getDeltaSphereRadius();
        		double kTension = be->kTension;
                U_i += _FFType.energy(kTension, d_radius);

        }
        else {

        		double d_radius = sb->getDeltaSphereRadius();
        		double kTension =be->kTension;
                U_i += _FFType.energy(kTension, d_radius, d);
            }
        }

        if(fabs(U_i) == numeric_limits<double>::infinity()
           || U_i != U_i || U_i < -1.0)
        {return -1;
        }
        else
            U += U_i;

        return U;

}



template <class BTensionType>
void BoundaryTension<BTensionType>::computeForces() {

    for (auto be: BoundaryElement::getBoundaryElements()) {
    		SphereBoundaryElement* sb = (SphereBoundaryElement*)be;
            double d_radius = sb->getDeltaSphereRadius();
    		double kTension = be->kTension;


            _FFType.forces(sb, kTension, d_radius);

    }
}


template <class BTensionType>
void BoundaryTension<BTensionType>::computeForcesAux() {

	for (auto be: BoundaryElement::getBoundaryElements()) {
		SphereBoundaryElement* sb = (SphereBoundaryElement*)be;

		double d_radius = sb->getDeltaSphereRadius();
		double kTension = be->kTension;



            _FFType.forcesAux(sb, kTension, d_radius);
        }
}

///Template specializations
template double BoundaryTension<BoundaryTensionHarmonic>::computeEnergy(double d);
template void BoundaryTension<BoundaryTensionHarmonic>::computeForces();
template void BoundaryTension<BoundaryTensionHarmonic>::computeForcesAux();
