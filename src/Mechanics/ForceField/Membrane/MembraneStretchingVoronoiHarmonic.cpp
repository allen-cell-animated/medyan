#include "MembraneStretchingVoronoiHarmonic.h"

#include "Bead.h"

#include "MathFunctions.h"

using namespace mathfunc;

// TODO: implement it

double MembraneStretchingHarmonic::energy(Bead* b, std::vector<Bead*> 
                                          double kElastic, double eqArea){
    // kElastic is the elastic modulus, which is independent of the actual eqArea

    double dist = currentArea - eqArea;
    
    return 0.5 * kElastic * dist * dist / eqArea;
    
}

double MembraneStretchingHarmonic::energy(const std::array<Bead*, 3>& b,
                                          double kElastic, double eqArea, double d){

    double distStretched = areaTriangleStretched(b[0]->coordinate, b[0]->force,
                                                 b[1]->coordinate, b[1]->force,
                                                 b[2]->coordinate, b[2]->force,
                                                 d) - eqArea;
    return 0.5 * kElastic * distStretched * distStretched / eqArea;
}

void MembraneStretchingHarmonic::forces(const std::array<Bead*, 3>& b,
                                        double kElastic, double eqArea ){
    
    double area = areaTriangle(b[0]->coordinate, b[1]->coordinate, b[2]->coordinate);
    double f0 = kElastic * (area - eqArea) / eqArea;

    // Actually area is calculated again in areaGradient function. There might be a way to optimize it...
    auto areaGradient = areaTriangleGradient(b[0]->coordinate, b[1]->coordinate, b[2]->coordinate);
    // force = kElastic * (area - eqArea) / eqArea * d_A

    for(int idxBead = 0; idxBead < 3; ++idxBead) {
        for(int idxCoord = 0; idxCoord < 3; ++idxCoord) {
            b[idxBead]->force[idxCoord] += f0 * areaGradient[idxBead][idxCoord];
        }
    }
    // TODO: test the results
}

void MembraneStretchingHarmonic::forcesAux(const std::array<Bead*, 3>& b,
                                           double kElastic, double eqArea ){
    
    double area = areaTriangle(b[0]->coordinate, b[1]->coordinate, b[2]->coordinate);
    double f0 = kElastic * (area - eqArea) / eqArea;

    auto areaGradient = areaTriangleGradient(b[0]->coordinate, b[1]->coordinate, b[2]->coordinate);

    for(int idxBead = 0; idxBead < 3; ++idxBead) {
        for(int idxCoord = 0; idxCoord < 3; ++idxCoord) {
            b[idxBead]->forceAux[idxCoord] += f0 * areaGradient[idxBead][idxCoord];
        }
    }
}

