/*
 * BoundaryTensionHarmonic.cpp
 *
 *  Created on: Jul 19, 2016
 *      Author: james
 */

#include "BoundaryElementImpl.h"
#include "BoundaryTensionHarmonic.h"


double BoundaryTensionHarmonic::energy(double kTension, double d_radius){

	return 0.5 * kTension * d_radius * d_radius;
}

double BoundaryTensionHarmonic::energy(double kTension, double d_radius, double d){

	return 0.5 * kTension * (d_radius + d) * (d_radius + d);

}

void BoundaryTensionHarmonic::forces(SphereBoundaryElement* b2, double kTension, double d_radius){

	double f0 = - kTension * d_radius;

	//initial total force on the boundary (boundary tension and exerted force by actins)
	b2->boundary_force +=  f0;

}

void BoundaryTensionHarmonic::forcesAux(SphereBoundaryElement* b2,double kTension, double d_radius){

	double f0 = - kTension * d_radius;

	//changing total force on the boundary
	b2->boundary_forceAux +=  f0;
}





