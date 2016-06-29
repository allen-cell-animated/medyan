
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v3.0
//
//  Copyright (2015)  Papoian Lab, University of Maryland
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
    
    double dist = twoPointDistance(b1->coordinate, b2->coordinate);
    
    double effd = dist - radius;
    
    double R = -effd / screenLength;
    return kRep * exp(R);
}

double BubbleCylinderRepulsionExp::energy(Bead* b1, Bead* b2, double radius,
                                          double kRep, double screenLength, double d) {
    
    double dist = twoPointDistanceStretched(b1->coordinate, b1->force,
                                            b2->coordinate, b2->force, d);
    double effd = dist - radius;
    
    double R = -effd / screenLength;
    return kRep * exp(R);
}

void BubbleCylinderRepulsionExp::forces(Bead* b1, Bead* b2, double radius,
                                        double kRep, double screenLength) {
    
    //get dist
    double dist = twoPointDistance(b1->coordinate, b2->coordinate);
    
    double effd = dist - radius;
    
    double R = -effd / screenLength;
    double f0 = kRep * exp(R) / screenLength;
    
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
    
    //get dist
    double dist = twoPointDistance(b1->coordinate, b2->coordinate);
    
    double effd = dist - radius;
    
    double R = -effd / screenLength;
    double f0 = kRep * exp(R) / screenLength;
    
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
    
    //get dist
    double dist = twoPointDistance(b1->coordinate, b2->coordinate);
    
    double effd = dist - radius;
    
    double R = -effd / screenLength;
    return kRep * exp(R) / screenLength;
}
