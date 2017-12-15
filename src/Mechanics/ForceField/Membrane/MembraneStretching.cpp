#include <array>

#include "MembraneStretching.h"

#include "Membrane.h"
#include "Triangle.h"
#include "Vertex.h"

#include "MembraneStretchingHarmonic.h"
#include "MembraneStretchingVoronoiHarmonic.h"

// Using the area of the Voronoi cells
template<>
double MembraneStretching<MembraneStretchingVoronoiHarmonic>::computeEnergy(double d) {
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
                // As a result, d is just a dummy variable due to its universality
                double areaStretched = v->getMVoronoiCell()->getStretchedArea();

                U_i += _FFType.energy(areaStretched, kElastic, eqArea, d);
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

template<>
void MembraneStretching<MembraneStretchingVoronoiHarmonic>::computeForces() {
    
    for (auto m: Membrane::getMembranes()) {
    
        for(Vertex* v : m->getVertexVector()){
            
            double kElastic = v->getMVoronoiCell()->getElasticModulus();
            double eqArea = v->getMVoronoiCell()->getEqArea();

            double area = v->getMVoronoiCell()->getArea();
            std::array<double, 3>& dArea = v->getMVoronoiCell()->getDArea();
            std::vector<std::array<double, 3>>& dNeighborArea = v->getMVoronoiCell()->getDNeighborArea();
           
            _FFType.forces(v, v->getNeighborVertices(), area, dArea, dNeighborArea, kElastic, eqArea);
        }
    }
}

template<>
void MembraneStretching<MembraneStretchingVoronoiHarmonic>::computeForcesAux() {
    
    for (auto m: Membrane::getMembranes()) {
    
        for(Vertex* v : m->getVertexVector()){
            
            double kElastic = v->getMVoronoiCell()->getElasticModulus();
            double eqArea = v->getMVoronoiCell()->getEqArea();

            double area = v->getMVoronoiCell()->getArea();
            std::array<double, 3>& dArea = v->getMVoronoiCell()->getDArea();
            std::vector<std::array<double, 3>>& dNeighborArea = v->getMVoronoiCell()->getDNeighborArea();
           
            _FFType.forcesAux(v, v->getNeighborVertices(), area, dArea, dNeighborArea, kElastic, eqArea);
        }
    }
}



// Using the areas of the triangles
template<>
double MembraneStretching<MembraneStretchingHarmonic>::computeEnergy(double d) {
    double U = 0;
    double U_i;

    for(auto m: Membrane::getMembranes()) {

        U_i = 0;

        if(d == 0.0){
            for(Triangle* it: m->getTriangleVector()){
                double kElastic = it->getMTriangle()->getElasticModulus();
                double eqArea = it->getMTriangle()->getEqArea();

                double area = it->getMTriangle()->getArea();

                // The calculation requires that the current area has already been calculated
                U_i += _FFType.energy(area, kElastic, eqArea);
            }

        } else {
            for(Triangle* it: m->getTriangleVector()){
                double kElastic = it->getMTriangle()->getElasticModulus();
                double eqArea = it->getMTriangle()->getEqArea();

                // The calculation requires that the current stretched area has been calculated
                double areaStretched = it->getMTriangle()->getStretchedArea();

                // Currently, d is a dummy variable, as the stretched areaStretched is already dependent on d.
                U_i += _FFType.energy(areaStretched, kElastic, eqArea, d);
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

template<>
void MembraneStretching<MembraneStretchingHarmonic>::computeForces() {
    
    for (auto m: Membrane::getMembranes()) {
    
        for(Triangle* it : m->getTriangleVector()){
            
            double kElastic = it->getMTriangle()->getElasticModulus();
            double eqArea = it->getMTriangle()->getEqArea();

            double area = it->getMTriangle()->getArea();
            std::array<std::array<double, 3>, 3>& dArea = it->getMTriangle()->getDArea();
           
            _FFType.forces(it->getVertices(), area, dArea, kElastic, eqArea);
        }
    }
}

template<>
void MembraneStretching<MembraneStretchingHarmonic>::computeForcesAux() {
    
    for (auto m: Membrane::getMembranes()) {
    
        for(Triangle* it : m->getTriangleVector()){
            
            double kElastic = it->getMTriangle()->getElasticModulus();
            double eqArea = it->getMTriangle()->getEqArea();

            double area = it->getMTriangle()->getArea();
            std::array<std::array<double, 3>, 3>& dArea = it->getMTriangle()->getDArea();
           
            _FFType.forcesAux(it->getVertices(), area, dArea, kElastic, eqArea);
        }
    }
}
