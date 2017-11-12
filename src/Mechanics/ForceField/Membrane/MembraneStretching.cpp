
#include "MembraneStretching.h"

#include "Membrane.h"
#include "Triangle.h"

template <class MembraneStretchingInteractionType>
double MembraneStretching<MembraneStretchingInteractionType>::computeEnergy(double d) {
    double U = 0;
    double U_i;

    for(auto m: Membrane::GetMembranes()) {

        U_i = 0;

        if(d == 0.0){
            for(Triangle* it: m->getTriangleVector()){
                double kElastic = it->getMTriangle()->getElasticModulus();
                double eqArea = it->getMTriangle()->getEqArea();

                // TODO: maybe I can calculate area first and store it in MTriangle
                U_i += _FFType.energy(it->getBeads(), kElastic, eqArea);
            }

        } else {
            for(Triangle* it: m->getTriangleVector()){
                double kElastic = it->getMTriangle()->getElasticModulus();
                double eqArea = it->getMTriangle()->getEqArea();

                U_i += _FFType.energy(it->getBeads(), kElastic, eqArea, d);
            }
        }

        if(fabs(U_i) == numeric_limits<double>::infinity()
           || U_i != U_i || U_i < -1.0) { // (U_i != U_i) => (U_i == NaN)
            
            //set culprit and return
            _membraneCulprit = m;
            
            return -1;
        }
        else
            U += U_i;

    }

    return U;
}

// TODO: implement others
