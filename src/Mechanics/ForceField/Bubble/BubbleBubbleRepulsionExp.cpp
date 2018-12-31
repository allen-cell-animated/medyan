
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v3.2.1
//
//  Copyright (2015-2018)  Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the MEDYAN web page for more information:
//  http://www.medyan.org
//------------------------------------------------------------------

#include "BubbleBubbleRepulsionExp.h"

#include "Bead.h"

#include "MathFunctions.h"

using namespace mathfunc;

double BubbleBubbleRepulsionExp::energy(Bead* b1, Bead* b2, double r1, double r2,
                                        double kRep, double screenLength, bool stretched) {

    const auto& c1 = stretched ? b1->getCoordinate<true>() : b1->getCoordinate<false>();
    const auto& c2 = stretched ? b2->getCoordinate<true>() : b2->getCoordinate<false>();

    double dist = twoPointDistance(c1, c2);
    
    double effd = dist - r1 - r2;
    
    double R = -effd / screenLength;
    return kRep * exp(R);
}

double BubbleBubbleRepulsionExp::energy(Bead* b1, Bead* b2, double r1, double r2,
                                        double kRep, double screenLength, double d) {
    
    double dist = twoPointDistanceStretched(b1->coordinate, b1->force,
                                            b2->coordinate, b2->force, d);
    double effd = dist - r1 - r2;
    
    double R = -effd / screenLength;
    return kRep * exp(R);
}

void BubbleBubbleRepulsionExp::forces(Bead* b1, Bead* b2, double r1, double r2,
                                      double kRep, double screenLength) {
    
    //get dist
    double dist = twoPointDistance(b1->coordinate, b2->coordinate);
    
    double effd = dist - r1 - r2;
    
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

void BubbleBubbleRepulsionExp::forcesAux(Bead* b1, Bead* b2, double r1, double r2,
                                         double kRep, double screenLength) {
    
    //get dist
    double dist = twoPointDistance(b1->coordinate, b2->coordinate);
    
    double effd = dist - r1 - r2;
    
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
