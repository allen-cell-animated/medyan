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

        if(d == 0.0) {
            const auto kElastic = m->getMMembrane()->getKElastic();
            const auto eqArea = m->getMMembrane()->getEqArea();

            double area = 0.0;
            for(const auto& v : m->getMesh().getVertices()) area += v.attr.gVertex.area;

            U_i = _FFType.energy(area, kElastic, eqArea); 

        } else {
            const auto kElastic = m->getMMembrane()->getKElastic();
            const auto eqArea = m->getMMembrane()->getEqArea();

            double sArea = 0.0;
            for(const auto& v : m->getMesh().getVertices()) sArea += v.attr.gVertex.sArea;

            U_i = _FFType.energy(sArea, kElastic, eqArea, d);
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

        const auto& mesh = m->getMesh();

        const auto kElastic = m->getMMembrane()->getKElastic();
        const auto eqArea = m->getMMembrane()->getEqArea();

        double area = 0.0;
        for(const auto& v : m->getMesh().getVertices()) area += v.attr.gVertex.area;

        const size_t numVertices = m->getMesh().getVertices().size();
        for(size_t vi = 0; vi < numVertices; ++vi) {
            const auto& v = m->getMesh().getVertices()[vi];
           
            _FFType.forces(v.attr.vertex, area, v.attr.gVertex.dArea, kElastic, eqArea);

            // Position of this vertex also affects neighbor vcell areas
            mesh.forEachHalfEdgeTargetingVertex(vi, [this, &mesh, &v, area, kElastic, eqArea](size_t hei) {
                const auto& dArea = mesh.getEdgeAttribute(hei).gHalfEdge.dNeighborArea;
                _FFType.forces(v.attr.vertex, area, dArea, kElastic, eqArea);
            });
        }
    }
}

template<>
void MembraneStretching<MembraneStretchingVoronoiHarmonic>::computeForcesAux() {
    
    for (auto m: Membrane::getMembranes()) {
    
        const auto& mesh = m->getMesh();

        const auto kElastic = m->getMMembrane()->getKElastic();
        const auto eqArea = m->getMMembrane()->getEqArea();

        double area = 0.0;
        for(const auto& v : m->getMesh().getVertices()) area += v.attr.gVertex.area;

        const size_t numVertices = m->getMesh().getVertices().size();
        for(size_t vi = 0; vi < numVertices; ++vi) {
            const auto& v = m->getMesh().getVertices()[vi];
           
            _FFType.forcesAux(v.attr.vertex, area, v.attr.gVertex.dArea, kElastic, eqArea);

            // Position of this vertex also affects neighbor vcell areas
            mesh.forEachHalfEdgeTargetingVertex(vi, [this, &mesh, &v, area, kElastic, eqArea](size_t hei) {
                const auto& dArea = mesh.getEdgeAttribute(hei).gHalfEdge.dNeighborArea;
                _FFType.forcesAux(v.attr.vertex, area, dArea, kElastic, eqArea);
            });
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

                double area = it->getGTriangle()->getArea();

                // The calculation requires that the current area has already been calculated
                U_i += _FFType.energy(area, kElastic, eqArea);
            }

        } else {
            for(Triangle* it: m->getTriangleVector()){
                double kElastic = it->getMTriangle()->getElasticModulus();
                double eqArea = it->getMTriangle()->getEqArea();

                // The calculation requires that the current stretched area has been calculated
                double areaStretched = it->getGTriangle()->getStretchedArea();

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

            double area = it->getGTriangle()->getArea();
            std::array<std::array<double, 3>, 3>& dArea = it->getGTriangle()->getDArea();
           
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

            double area = it->getGTriangle()->getArea();
            std::array<std::array<double, 3>, 3>& dArea = it->getGTriangle()->getDArea();
           
            _FFType.forcesAux(it->getVertices(), area, dArea, kElastic, eqArea);
        }
    }
}
