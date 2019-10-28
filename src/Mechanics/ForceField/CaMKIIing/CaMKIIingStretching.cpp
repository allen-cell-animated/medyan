
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

#include "CaMKIIingStretching.h"

#include "CaMKIIingStretchingHarmonic.h"

#include "CaMKIIingPoint.h"
#include "Cylinder.h"
#include "Bead.h"

#include "MathFunctions.h"
using namespace mathfunc;

template <class CStretchingInteractionType>
double CaMKIIingStretching<CStretchingInteractionType>::computeEnergy(double d) {
    
    double U = 0;
    double U_i;

    for (auto point : CaMKIIingPoint::getCaMKIIingPoints()) {
        for (auto bond : point->getBonds()) {

            Cylinder *b = (Cylinder *) get<0>(bond);

            Bead *b1 = b->getFirstBead();
            Bead *b2 = b->getSecondBead();
            Bead *b3 = point->getCaMKIICylinder()->getFirstBead();

            double kStretch = point->getMCaMKIIingPoint()->getStretchingConstant();
            double eqLength = point->getMCaMKIIingPoint()->getEqLength();
            double position = point->getPosition();
            if (d == 0.0)
                U_i = _FFType.energy(b1, b2, b3, position, kStretch, eqLength);
            else
                U_i = _FFType.energy(b1, b2, b3, position, kStretch, eqLength, d);

            if (fabs(U_i) == numeric_limits<double>::infinity() || U_i != U_i || U_i < -1.0) {
                //set culprit and return
                _camkiiingCulprit = point;
                return -1;
            } else
                U += U_i;
        }
    }
    return U;
}

template <class CStretchingInteractionType>
void CaMKIIingStretching<CStretchingInteractionType>::computeForces() {
	auto i = 0;

	for (auto point : CaMKIIingPoint::getCaMKIIingPoints()) {
		for (auto bond : point->getBonds()) {

			Cylinder *b = (Cylinder *) get<0>(bond);

			Bead *b1 = b->getFirstBead();
			Bead *b2 = b->getSecondBead();
			Bead *b3 = point->getCaMKIICylinder()->getFirstBead();

			double kStretch = point->getMCaMKIIingPoint()->getStretchingConstant();
			double eqLength = point->getMCaMKIIingPoint()->getEqLength();
			double position = point->getPosition();

			_FFType.forces(b1, b2, b3, position, kStretch, eqLength);
		}
	}
}


template <class CStretchingInteractionType>
void CaMKIIingStretching<CStretchingInteractionType>::computeForcesAux() {
    
    for (auto point : CaMKIIingPoint::getCaMKIIingPoints()) {
        for(auto bond : point->getBonds()) {

            Cylinder* b = (Cylinder*) get<0>(bond);

            Bead *b1 = b->getFirstBead();
            Bead *b2 = b->getSecondBead();
            Bead *b3 = point->getCaMKIICylinder()->getFirstBead();

            double kStretch = point->getMCaMKIIingPoint()->getStretchingConstant();
            double eqLength = point->getMCaMKIIingPoint()->getEqLength();
            double position = point->getPosition();
            _FFType.forcesAux(b1, b2, b3, position, kStretch, eqLength);
        }
    }
}


///Template specializations
template double CaMKIIingStretching<CaMKIIingStretchingHarmonic>::computeEnergy(double d);
template void CaMKIIingStretching<CaMKIIingStretchingHarmonic>::computeForces();
template void CaMKIIingStretching<CaMKIIingStretchingHarmonic>::computeForcesAux();
