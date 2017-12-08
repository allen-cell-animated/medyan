
#include "MembraneStretching.h"

#include "Membrane.h"
#include "Triangle.h"
#include "Vertex.h"

#include "MembraneStretchingHarmonic.h"
#include "MembraneStretchingVoronoiHarmonic.h"

// Using the area of the Voronoi cells
template double MembraneStretching<MembraneStretchingVoronoiHarmonic>::computeEnergy(double d) {
    double U = 0;
    double U_i;

    for(auto m: Membrane::getMembranes()) {
        U_i = 0;

        if(d == 0.0) {
            for(Vertex* v: m->getVertexVector()) {
                double kElastic = v->getMVoronoiCell()->getElasticModulus();
                double eqArea = v->getMVoronoiCell()->getEqArea();

                // The calculation requires that the current area has already been calculated
                double area = v->getMVoronoiCell()->getArea();

                U_i += _FFType.energy(area, kElastic, eqArea); 
            }
        } else {
            for(Vertex *v: m->getVertexVector()) {
                double kElastic = v->getMVoronoiCell()->getElasticModulus();
                double eqArea = v->getMVoronoiCell()->getEqArea();

                // The calculation requires that the current stretched area has already been calculated
                double area = v->getMVoronoiCell()->getTempArea(); // TODO: implement this function

                U_i += _FFType.energy(area, kElastic, eqArea, d);
            }
        }

        if(fabs(U_i) == numeric_limits<double>::infinity()
            || U_i != U_i || U_i < -1.0) {
            _membraneCulprit = m;
            return -1;
        } else
            U += U_i;
        
    }

    return U;
}

// Using the areas of the triangles
template double MembraneStretching<MembraneStretchingHarmonic>::computeEnergy(double d) {
    double U = 0;
    double U_i;

    for(auto m: Membrane::getMembranes()) {

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

template <class MembraneStretchingInteractionType>
void MembraneStretching<MembraneStretchingInteractionType>::computeForces() {
    
    for (auto m: Membrane::getMembranes()) {
    
        for(Triangle* it : f->getTriangleVector()){
            
            double kElastic =it->getMTriangle()->getElasticModulus();
            double eqArea = it->getMTriangle()->getEqArea();
           
            _FFType.forces(it->getBeads(), kElastic, eqArea);
        }
    }
}

template <class MembraneStretchingInteractionType>
void MembraneStretching<MembraneStretchingInteractionType>::computeForcesAux() {
    
    for (auto m: Membrane::getMembranes()) {
    
        for(Triangle* it : f->getTriangleVector()){
            
            double kElastic =it->getMTriangle()->getElasticModulus();
            double eqArea = it->getMTriangle()->getEqArea();
           
            _FFType.forcesAux(it->getBeads(), kElastic, eqArea);
        }
    }
}

// Template specializations
template double MembraneStretching<MembraneStretchingHarmonic>::computeEnergy(double d);
template void MembraneStretching<MembraneStretchingHarmonic>::computeForces();
template void MembraneStretching<MembraneStretchingHarmonic>::computeForcesAux();
