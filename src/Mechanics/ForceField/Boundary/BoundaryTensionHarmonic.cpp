/*
 * BoundaryTensionHarmonic.cpp
 *
 *  Created on: Jul 19, 2016
 *      Author: james
 */

#include "BoundaryElementImpl.h"
#include "BoundaryTensionHarmonic.h"


double BoundaryTensionHarmonic::energy(double kTension, double radius){

	return 0.5 * kTension * radius * radius;
}

double BoundaryTensionHarmonic::energy(double kTension, double radius, double d){

	return 0.5 * kTension * (radius + d) * (radius + d);

}

void BoundaryTensionHarmonic::forces(SphereBoundaryElement* b2, double kTension, double radius){

	double f0 = kTension * radius;

	//initial boundary tension
	b2->boundary_tension +=  f0;

}

void BoundaryTensionHarmonic::forcesAux(SphereBoundaryElement* b2,double kTension, double radius){

	double f0 = kTension * radius;

	//changing boundary tension
	b2->boundary_tensionaux +=  f0;
}





