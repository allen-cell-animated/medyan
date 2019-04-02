
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

template <class BStretchingInteractionType>
double CaMKIIingStretching<BStretchingInteractionType>::computeEnergy(double d) {
    
    double U = 0;
    double U_i;
    
    for (auto b: CaMKIIingPoint::getCaMKIIingPoints()) {
        
        Bead* b1 = b->getFirstCylinder()->getFirstBead();
        Bead* b2 = b->getFirstCylinder()->getSecondBead();
        Bead* b3 = b->getSecondCylinder()->getFirstBead();
        
        double kStretch = b->getMCaMKIIingPoint()->getStretchingConstant();
        double eqLength = b->getMCaMKIIingPoint()->getEqLength();
        double position = b->getPosition();
// TODO fix to iterate over each of the _bonds according to coordination number
#if 0
        if (d == 0.0)
            U_i = _FFType.energy(b1, b2, b3, position, kStretch, eqLength);
        else
            U_i = _FFType.energy(b1, b2, b3, position, kStretch, eqLength, d);
        
        if(fabs(U_i) == numeric_limits<double>::infinity()
           || U_i != U_i || U_i < -1.0) {
            
            //set culprit and return
            _camkiiingCulprit = b;
            
            return -1;
        }
        else
            U += U_i;
#endif
    }
    return U;
}

template <class BStretchingInteractionType>
void CaMKIIingStretching<BStretchingInteractionType>::computeForces() {
    auto i=0;

    cout << "CAMKII" << __LINE__ << endl;
    for (auto b: CaMKIIingPoint::getCaMKIIingPoints()) {
        i++;
        cout << "CAMKII" << i << " " << __LINE__ << endl;
        Bead* b1 = b->getFirstCylinder()->getFirstBead();
        cout << "CAMKII" << i << " " << __LINE__ << endl;
        Bead* b2 = b->getFirstCylinder()->getSecondBead();
        cout << "CAMKII" << i << " " << __LINE__ << endl;
        Bead* b3 = b->getSecondCylinder()->getFirstBead();
        cout << "CAMKII" << i << " " << __LINE__ << endl;
        Bead* b4 = b->getSecondCylinder()->getSecondBead();

        std::cout<<i<<" "<<b->getFirstCylinder()->getID()<<" "<<twoPointDistance(b1->coordinate, b2->coordinate)<<" "<<b->getSecondCylinder()->getID()<<" "<<twoPointDistance(b3->coordinate, b4->coordinate)<<endl;

        double kStretch = b->getMCaMKIIingPoint()->getStretchingConstant();
        cout << "CAMKII" << b << " " << __LINE__ << endl;
        double eqLength = b->getMCaMKIIingPoint()->getEqLength();
        cout << "CAMKII" << b << " " << __LINE__ << endl;
        double position = b->getPosition();
        // TODO fix to iterate over each of the _bonds according to coordination number
        //_FFType.forces(b1, b2, b3, position, kStretch, eqLength);
    }
}


template <class BStretchingInteractionType>
void CaMKIIingStretching<BStretchingInteractionType>::computeForcesAux() {
    
    for (auto b: CaMKIIingPoint::getCaMKIIingPoints()) {
    	cout << "CAMKII" << b << " " << __LINE__ << endl;
        Bead* b1 = b->getFirstCylinder()->getFirstBead();
        cout << "CAMKII" << b << " " << __LINE__ << endl;
        Bead* b2 = b->getFirstCylinder()->getSecondBead();
        cout << "CAMKII" << b << " " << __LINE__ << endl;
        Bead* b3 = b->getSecondCylinder()->getFirstBead();
        cout << "CAMKII" << b << " " << __LINE__ << endl;
        
        double kStretch = b->getMCaMKIIingPoint()->getStretchingConstant();
        cout << "CAMKII" << b << " " << __LINE__ << endl;
        double eqLength = b->getMCaMKIIingPoint()->getEqLength();
        cout << "CAMKII" << b << " " << __LINE__ << endl;
        double position = b->getPosition();
        // TODO fix to iterate over each of the _bonds according to coordination number
        //_FFType.forcesAux(b1, b2, b3, position, kStretch, eqLength);
    }
}


///Template specializations
template double
CaMKIIingStretching<CaMKIIingStretchingHarmonic>::computeEnergy(double d);
template void CaMKIIingStretching<CaMKIIingStretchingHarmonic>::computeForces();
template void CaMKIIingStretching<CaMKIIingStretchingHarmonic>::computeForcesAux();
