#include "MembraneStretchingHarmonic.h"

#include "Bead.h"

#include "MathFunctions.h"

using namespace mathfunc;

double MembraneStretchingHarmonic::energy(const std::array<Bead*, 3>& b,
                                          double kElastic, double eqArea){
    // kElastic is the elastic modulus, which is independent of the actual eqArea

    double dist = areaTriangle(b[0]->coordinate, b[1]->coordinate, b[2]->coordinate) - eqArea;
    
    return 0.5 * kStretch* dist * dist / eqArea;
    
}

double MembraneStretchingHarmonic::energy(Bead* b1, Bead* b2,
                                          double kStretch, double eqLength, double d){

    double distStretched = twoPointDistanceStretched(b1->coordinate,
                                                     b1->force, b2->coordinate,
                                                     b2->force, d) - eqLength;
    return 0.5 * kStretch * distStretched * distStretched;
}

void MembraneStretchingHarmonic::forces(Bead* b1, Bead* b2,
                                        double kStretch, double eqLength ){
    
    double dist = twoPointDistance( b1->coordinate, b2->coordinate);
    double invL = 1 / dist;
    
    double f0 = kStretch * ( dist - eqLength ) * invL;
    
    //force on i
    b2->force[0] +=  f0 * ( b1->coordinate[0] - b2->coordinate[0] );
    b2->force[1] +=  f0 * ( b1->coordinate[1] - b2->coordinate[1] );
    b2->force[2] +=  f0 * ( b1->coordinate[2] - b2->coordinate[2] );
    
    // force i-1
    b1->force[0] +=  f0 * ( b2->coordinate[0] - b1->coordinate[0] );
    b1->force[1] +=  f0 * ( b2->coordinate[1] - b1->coordinate[1] );
    b1->force[2] +=  f0 * ( b2->coordinate[2] - b1->coordinate[2] );
}

void MembraneStretchingHarmonic::forcesAux(Bead* b1, Bead* b2,
                                           double kStretch, double eqLength ){
    
    double dist = twoPointDistance( b1->coordinate, b2->coordinate);
    double invL = 1 / dist;
    double f0 = kStretch * ( dist - eqLength ) * invL;

    //force on i
    b2->forceAux[0] +=  f0 * ( b1->coordinate[0] - b2->coordinate[0] );
    b2->forceAux[1] +=  f0 * ( b1->coordinate[1] - b2->coordinate[1] );
    b2->forceAux[2] +=  f0 * ( b1->coordinate[2] - b2->coordinate[2] );
    
    // force i-1
    b1->forceAux[0] +=  f0 * ( b2->coordinate[0] - b1->coordinate[0] );
    b1->forceAux[1] +=  f0 * ( b2->coordinate[1] - b1->coordinate[1] );
    b1->forceAux[2] +=  f0 * ( b2->coordinate[2] - b1->coordinate[2] );
}

