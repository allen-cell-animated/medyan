
#include "MembraneStretching.h"

#include "Membrane.h"
#include "Triangle.h"

template <class MembraneStretchingInteractionType>
double MembraneStretching<MembraneStretchingInteractionType>::computeEnergy(double d) {
    double U = 0;
    double U_i;

    for(auto m: Membrane::GetMembranes()) {

        U_i = 0;

        // TODO: content
        if(d == 0.0){

        }
    }

    return U;
}
