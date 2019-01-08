
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

#include <cmath>
#include <math.h>

#include "BranchingPositionCosine.h"

#include "Bead.h"

#include "MathFunctions.h"


using namespace mathfunc;

double BranchingPositionCosine::energy(Bead* b1, Bead* b2, Bead* b3,
                                       double kPosition, double position, bool stretched){

    const auto& c1 = stretched ? b1->getCoordinate<true>() : b1->getCoordinate<false>();
    const auto& c2 = stretched ? b2->getCoordinate<true>() : b2->getCoordinate<false>();
    const auto& c3 = stretched ? b3->getCoordinate<true>() : b3->getCoordinate<false>();

    double X = sqrt(scalarProduct(
    midPointCoordinate(c1, c2, position),c2,
    midPointCoordinate(c1, c2, position),c2));
    
    double D = sqrt(scalarProduct(
    midPointCoordinate(c1, c2, position), c3,
    midPointCoordinate(c1, c2, position), c3));
    
    double XD = X*D;
    
    double xd = scalarProduct(
    midPointCoordinate(c1, c2, position),c2,
    midPointCoordinate(c1, c2, position), c3);
    
    double theta = safeacos(xd / XD);
    double eqTheta = 0.5*M_PI;
    double dtheta = theta-eqTheta;
    
    double U = kPosition * ( 1 - cos(dtheta) );
    
    return U;
    
}

void BranchingPositionCosine::forces(Bead* b1, Bead* b2, Bead* b3,
                                     double kPosition, double position ){
    
    double X = sqrt(scalarProduct(
    midPointCoordinate(b1->coordinate, b2->coordinate, position),b2->coordinate,
    midPointCoordinate(b1->coordinate, b2->coordinate, position),b2->coordinate));
    
    double D = sqrt(scalarProduct(
    midPointCoordinate(b1->coordinate, b2->coordinate, position), b3->coordinate,
    midPointCoordinate(b1->coordinate, b2->coordinate, position), b3->coordinate));
    
    double xd = scalarProduct(
    midPointCoordinate(b1->coordinate, b2->coordinate, position),b2->coordinate,
    midPointCoordinate(b1->coordinate, b2->coordinate, position), b3->coordinate);
    
    //invL = 1/L;
    double invX = 1/X;
    double invD = 1/D;
    double A = invX*invD;
    double B = invX*invX;
    double C = invD*invD;
    
    double theta = safeacos(xd * A);
    double eqTheta = 0.5*M_PI;
    double dtheta = theta-eqTheta;
    
    double k =  kPosition * A * sin(dtheta)/sin(theta);
    
    //
    
    b1->force[0] +=  k * (1-position)* (- (1-position)*(b2->coordinate[0] - b1->coordinate[0]) - (b3->coordinate[0] - (1-position)*b1->coordinate[0] - position*b2->coordinate[0]) + xd *(B*(1-position)*(b2->coordinate[0] - b1->coordinate[0]) + C*(b3->coordinate[0] - (1-position)*b1->coordinate[0] - position*b2->coordinate[0])) );
    
    b1->force[1] +=  k * (1-position)* (- (1-position)*(b2->coordinate[1] - b1->coordinate[1]) - (b3->coordinate[1] - (1-position)*b1->coordinate[1] - position*b2->coordinate[1]) + xd *(B*(1-position)*(b2->coordinate[1] - b1->coordinate[1]) + C*(b3->coordinate[1] - (1-position)*b1->coordinate[1] - position*b2->coordinate[1])) );
    
    b1->force[2] +=  k * (1-position)* (- (1-position)*(b2->coordinate[2] - b1->coordinate[2]) - (b3->coordinate[2] - (1-position)*b1->coordinate[2] - position*b2->coordinate[2]) + xd *(B*(1-position)*(b2->coordinate[2] - b1->coordinate[2]) + C*(b3->coordinate[2] - (1-position)*b1->coordinate[2] - position*b2->coordinate[2])) );
    
    //
    
    b2->force[0] +=  k * (- position*(1-position)*(b2->coordinate[0] - b1->coordinate[0]) + (1-position)*(b3->coordinate[0]- (1-position)*b1->coordinate[0] - position*b2->coordinate[0]) + xd *( (1-position)*B*(1-position)*(b2->coordinate[0] - b1->coordinate[0]) - position*C*(b3->coordinate[0] - (1-position)*b1->coordinate[0] - position*b2->coordinate[0])) );
    
    b2->force[1] +=  k * (- position*(1-position)*(b2->coordinate[1] - b1->coordinate[1]) + (1-position)*(b3->coordinate[1]- (1-position)*b1->coordinate[1] - position*b2->coordinate[1]) + xd *( (1-position)*B*(1-position)*(b2->coordinate[1] - b1->coordinate[1]) - position*C*(b3->coordinate[1] - (1-position)*b1->coordinate[1] - position*b2->coordinate[1])) );
    
    b2->force[2] +=  k * (- position*(1-position)*(b2->coordinate[2] - b1->coordinate[2]) + (1-position)*(b3->coordinate[2]- (1-position)*b1->coordinate[2] - position*b2->coordinate[2]) + xd *( (1-position)*B*(1-position)*(b2->coordinate[2] - b1->coordinate[2]) - position*C*(b3->coordinate[2] - (1-position)*b1->coordinate[2] - position*b2->coordinate[2])) );
    
    //
    
    b3->force[0] +=  k * ( (1-position)*(b2->coordinate[0] - b1->coordinate[0]) - xd * C*(b3->coordinate[0] - (1-position)*b1->coordinate[0] - position*b2->coordinate[0]) );
    b3->force[1] +=  k * ( (1-position)*(b2->coordinate[1] - b1->coordinate[1]) - xd * C*(b3->coordinate[1] - (1-position)*b1->coordinate[1] - position*b2->coordinate[1]) );
    b3->force[2] +=  k * ( (1-position)*(b2->coordinate[2] - b1->coordinate[2]) - xd * C*(b3->coordinate[2] - (1-position)*b1->coordinate[2] - position*b2->coordinate[2]) );

    
}

void BranchingPositionCosine::forcesAux(Bead* b1, Bead* b2, Bead* b3,
                                     double kPosition, double position ){
    
    double X = sqrt(scalarProduct(
    midPointCoordinate(b1->coordinate, b2->coordinate, position),b2->coordinate,
    midPointCoordinate(b1->coordinate, b2->coordinate, position),b2->coordinate));
    
    double D = sqrt(scalarProduct(
    midPointCoordinate(b1->coordinate, b2->coordinate, position), b3->coordinate,
    midPointCoordinate(b1->coordinate, b2->coordinate, position), b3->coordinate));
    
    double xd = scalarProduct(
    midPointCoordinate(b1->coordinate, b2->coordinate, position),b2->coordinate,
    midPointCoordinate(b1->coordinate, b2->coordinate, position), b3->coordinate);
    
    //invL = 1/L;
    double invX = 1/X;
    double invD = 1/D;
    double A = invX*invD;
    double B = invX*invX;
    double C = invD*invD;
    
    double theta = safeacos(xd * A);
    double eqTheta = 0.5*M_PI;
    double dtheta = theta-eqTheta;
    
    double k =  kPosition * A * sin(dtheta)/sin(theta);
    
    //
    
    b1->forceAux[0] +=  k * (1-position)* (- (1-position)*(b2->coordinate[0] - b1->coordinate[0]) - (b3->coordinate[0] - (1-position)*b1->coordinate[0] - position*b2->coordinate[0]) + xd *(B*(1-position)*(b2->coordinate[0] - b1->coordinate[0]) + C*(b3->coordinate[0] - (1-position)*b1->coordinate[0] - position*b2->coordinate[0])) );
    
    b1->forceAux[1] +=  k * (1-position)* (- (1-position)*(b2->coordinate[1] - b1->coordinate[1]) - (b3->coordinate[1] - (1-position)*b1->coordinate[1] - position*b2->coordinate[1]) + xd *(B*(1-position)*(b2->coordinate[1] - b1->coordinate[1]) + C*(b3->coordinate[1] - (1-position)*b1->coordinate[1] - position*b2->coordinate[1])) );
    
    b1->forceAux[2] +=  k * (1-position)* (- (1-position)*(b2->coordinate[2] - b1->coordinate[2]) - (b3->coordinate[2] - (1-position)*b1->coordinate[2] - position*b2->coordinate[2]) + xd *(B*(1-position)*(b2->coordinate[2] - b1->coordinate[2]) + C*(b3->coordinate[2] - (1-position)*b1->coordinate[2] - position*b2->coordinate[2])) );
    
    //
    
    b2->forceAux[0] +=  k * (- position*(1-position)*(b2->coordinate[0] - b1->coordinate[0]) + (1-position)*(b3->coordinate[0]- (1-position)*b1->coordinate[0] - position*b2->coordinate[0]) + xd *( (1-position)*B*(1-position)*(b2->coordinate[0] - b1->coordinate[0]) - position*C*(b3->coordinate[0] - (1-position)*b1->coordinate[0] - position*b2->coordinate[0])) );
    
    b2->forceAux[1] +=  k * (- position*(1-position)*(b2->coordinate[1] - b1->coordinate[1]) + (1-position)*(b3->coordinate[1]- (1-position)*b1->coordinate[1] - position*b2->coordinate[1]) + xd *( (1-position)*B*(1-position)*(b2->coordinate[1] - b1->coordinate[1]) - position*C*(b3->coordinate[1] - (1-position)*b1->coordinate[1] - position*b2->coordinate[1])) );
    
    b2->forceAux[2] +=  k * (- position*(1-position)*(b2->coordinate[2] - b1->coordinate[2]) + (1-position)*(b3->coordinate[2]- (1-position)*b1->coordinate[2] - position*b2->coordinate[2]) + xd *( (1-position)*B*(1-position)*(b2->coordinate[2] - b1->coordinate[2]) - position*C*(b3->coordinate[2] - (1-position)*b1->coordinate[2] - position*b2->coordinate[2])) );
    
    //
    
    b3->forceAux[0] +=  k * ( (1-position)*(b2->coordinate[0] - b1->coordinate[0]) - xd * C*(b3->coordinate[0] - (1-position)*b1->coordinate[0] - position*b2->coordinate[0]) );
    b3->forceAux[1] +=  k * ( (1-position)*(b2->coordinate[1] - b1->coordinate[1]) - xd * C*(b3->coordinate[1] - (1-position)*b1->coordinate[1] - position*b2->coordinate[1]) );
    b3->forceAux[2] +=  k * ( (1-position)*(b2->coordinate[2] - b1->coordinate[2]) - xd * C*(b3->coordinate[2] - (1-position)*b1->coordinate[2] - position*b2->coordinate[2]) );
    
    
}

