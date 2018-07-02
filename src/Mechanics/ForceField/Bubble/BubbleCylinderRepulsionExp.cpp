
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

#include "BubbleCylinderRepulsionExp.h"

#include "Bead.h"

#include "MathFunctions.h"

using namespace mathfunc;

double BubbleCylinderRepulsionExp::energy(Bead* b1, Bead* b2, double radius,
                                          double kRep, double screenLength) {
    
    if (tau()<=60.5) {
    kRep=40;
	} 
    else if (tau()<=120.5) {
    kRep=80;
	}
    else if (tau()<=180.5) {
    kRep=120;
	}
    else if (tau()<=240.5) {
    kRep=200;
	}
    else if (tau()<=300.5) {
    kRep=280;
	}
    else {
    kRep=360;
	}
    
    double dist = twoPointDistance(b1->coordinate, b2->coordinate);
    
    double effd = dist - radius;
    
    double R = -effd / screenLength;

    return 361-kRep * exp(R);
    
}

double BubbleCylinderRepulsionExp::energy(Bead* b1, Bead* b2, double radius,
                                          double kRep, double screenLength, double d) {
    
    if (tau()<=60.5) {
    kRep=40;
	} 
    else if (tau()<=120.5) {
    kRep=80;
	}
    else if (tau()<=180.5) {
    kRep=120;
	}
    else if (tau()<=240.5) {
    kRep=200;
	}
    else if (tau()<=300.5) {
    kRep=280;
	}
    else {
    kRep=360;
	}
    
    double dist = twoPointDistanceStretched(b1->coordinate, b1->force,
                                            b2->coordinate, b2->force, d);
    double effd = dist - radius;
    
    double R = -effd / screenLength;

    return 361-kRep * exp(R);

}

void BubbleCylinderRepulsionExp::forces(Bead* b1, Bead* b2, double radius,
                                        double kRep, double screenLength) {
    
    if (tau()<=60.5) {
    kRep=40;
	} 
    else if (tau()<=120.5) {
    kRep=80;
	}
    else if (tau()<=180.5) {
    kRep=120;
	}
    else if (tau()<=240.5) {
    kRep=200;
	}
    else if (tau()<=300.5) {
    kRep=280;
	}
    else {
    kRep=360;
	}
    
    //get dist
    double dist = twoPointDistance(b1->coordinate, b2->coordinate);
    
    double effd = dist - radius;
    
    double R = -effd / screenLength;

    double f0 = -kRep * exp(R) / screenLength;
            
    //get norm
    auto norm = normalizedVector(twoPointDirection(b1->coordinate, b2->coordinate));
    
    b1->force[0] += - f0 *norm[0];
    b1->force[1] += - f0 *norm[1];
    b1->force[2] += - f0 *norm[2];
    
    b2->force[0] += f0 *norm[0];
    b2->force[1] += f0 *norm[1];
    b2->force[2] += f0 *norm[2];
}

void BubbleCylinderRepulsionExp::forcesAux(Bead* b1, Bead* b2, double radius,
                                           double kRep, double screenLength) {
    if (tau()<=60.5) {
    kRep=40;
	} 
    else if (tau()<=120.5) {
    kRep=80;
	}
    else if (tau()<=180.5) {
    kRep=120;
	}
    else if (tau()<=240.5) {
    kRep=200;
	}
    else if (tau()<=300.5) {
    kRep=280;
	}
    else {
    kRep=360;
	}
    
    //get dist
    double dist = twoPointDistance(b1->coordinate, b2->coordinate);
    
    double effd = dist - radius;
    
    double R = -effd / screenLength;

    double f0 = -kRep * exp(R) / screenLength;

    //get norm
    auto norm = normalizedVector(twoPointDirection(b1->coordinate, b2->coordinate));
    
    b1->force[0] += - f0 *norm[0];
    b1->force[1] += - f0 *norm[1];
    b1->force[2] += - f0 *norm[2];
    
    b2->force[0] += f0 *norm[0];
    b2->force[1] += f0 *norm[1];
    b2->force[2] += f0 *norm[2];
    
}

double BubbleCylinderRepulsionExp::loadForces(Bead* b1, Bead* b2, double radius,
                                              double kRep, double screenLength) {
    
    if (tau()<=60.5) {
    kRep=40;
	} 
    else if (tau()<=120.5) {
    kRep=80;
	}
    else if (tau()<=180.5) {
    kRep=120;
	}
    else if (tau()<=240.5) {
    kRep=200;
	}
    else if (tau()<=300.5) {
    kRep=280;
	}
    else {
    kRep=360;
	}
    
    //get dist
    double dist = twoPointDistance(b1->coordinate, b2->coordinate);
    
    double effd = dist - radius;
    
    double R = -effd / screenLength;

    return -kRep * exp(R) / screenLength;

}
