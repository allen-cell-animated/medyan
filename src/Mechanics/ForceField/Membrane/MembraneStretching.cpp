#include <array>

#include "MembraneStretching.h"

#include "Membrane.hpp"
#include "Triangle.h"
#include "Vertex.h"

#include "MembraneStretchingHarmonic.h"
#include "MembraneStretchingVoronoiHarmonic.h"

// Using the area of the Voronoi cells
template<>
double MembraneStretching<MembraneStretchingVoronoiHarmonic>::computeEnergy(bool stretched) {
    double U = 0;
    double U_i;

    for(auto m: Membrane::getMembranes()) {

        const auto kElastic = m->getMMembrane()->getKElastic();
        const auto eqArea = m->getMMembrane()->getEqArea();

        double area = 0.0;
        for(const auto& v : m->getMesh().getVertices())
            area += stretched ? v.attr.gVertex.sArea : v.attr.gVertex.area;

        U_i = _FFType.energy(area, kElastic, eqArea); 

        if(fabs(U_i) == numeric_limits<double>::infinity()
            || U_i != U_i || U_i < -1.0) {
            _membraneCulprit = m;
            return -1;
        } else
            U += U_i;
        
    }

    return U;
}

// Force on each vertex is calculated one-time using derivative of vcell area and
// the derivative of neighbor vcell areas on the center vertex.
template<>
void MembraneStretching<MembraneStretchingVoronoiHarmonic>::computeForces() {
    
    for (auto m: Membrane::getMembranes()) {

        const auto& mesh = m->getMesh();

        const auto kElastic = m->getMMembrane()->getKElastic();
        const auto eqArea = m->getMMembrane()->getEqArea();

        double area = 0.0;
        for(const auto& v : mesh.getVertices()) area += v.attr.gVertex.area;

        const size_t numVertices = mesh.getVertices().size();
        for(size_t vi = 0; vi < numVertices; ++vi) {
            const auto& v = mesh.getVertices()[vi];
           
            _FFType.forces(v.attr.vertex, area, v.attr.gVertex.dArea, kElastic, eqArea);

            // Position of this vertex also affects neighbor vcell areas
            mesh.forEachHalfEdgeTargetingVertex(vi, [this, &mesh, &v, area, kElastic, eqArea](size_t hei) {
                const auto& dArea = mesh.getHalfEdgeAttribute(hei).gHalfEdge.dNeighborArea;
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
        for(const auto& v : mesh.getVertices()) area += v.attr.gVertex.area;

        const size_t numVertices = mesh.getVertices().size();
        for(size_t vi = 0; vi < numVertices; ++vi) {
            const auto& v = mesh.getVertices()[vi];
           
            _FFType.forcesAux(v.attr.vertex, area, v.attr.gVertex.dArea, kElastic, eqArea);

            // Position of this vertex also affects neighbor vcell areas
            mesh.forEachHalfEdgeTargetingVertex(vi, [this, &mesh, &v, area, kElastic, eqArea](size_t hei) {
                const auto& dArea = mesh.getHalfEdgeAttribute(hei).gHalfEdge.dNeighborArea;
                _FFType.forcesAux(v.attr.vertex, area, dArea, kElastic, eqArea);
            });
        }
    }
}



// Using the areas of the triangles
template<>
double MembraneStretching<MembraneStretchingHarmonic>::computeEnergy(bool stretched) {
    double U = 0;
    double U_i;

    for(auto m: Membrane::getMembranes()) {

        const auto kElastic = m->getMMembrane()->getKElastic();
        const auto eqArea = m->getMMembrane()->getEqArea();

        double area = 0.0;
        for(const auto& t : m->getMesh().getTriangles())
            area += stretched ? t.attr.gTriangle.sArea : t.attr.gTriangle.area;

        U_i = _FFType.energy(area, kElastic, eqArea); 

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

// Currently force calculation using triangles are different with the one using vcells.
// Using triangles, looping through triangles and forces are accumulated on the vertices.
template<>
void MembraneStretching<MembraneStretchingHarmonic>::computeForces() {
    
    for (auto m: Membrane::getMembranes()) {
    
        const auto& mesh = m->getMesh();

        const auto kElastic = m->getMMembrane()->getKElastic();
        const auto eqArea = m->getMMembrane()->getEqArea();

        double area = 0.0;
        for(const auto& t : mesh.getTriangles()) area += t.attr.gTriangle.area;

        const size_t numTriangles = mesh.getTriangles().size();
        for(size_t ti = 0; ti < numTriangles; ++ti) {
            mesh.forEachHalfEdgeInTriangle(ti, [this, &mesh, area, kElastic, eqArea](size_t hei) {
                const auto& dArea = mesh.getHalfEdgeAttribute(hei).gHalfEdge.dTriangleArea;
                _FFType.forces(mesh.getVertexAttribute(mesh.target(hei)).vertex, area, dArea, kElastic, eqArea);
            });
        }
    }
}

template<>
void MembraneStretching<MembraneStretchingHarmonic>::computeForcesAux() {
    
    for (auto m: Membrane::getMembranes()) {
    
        const auto& mesh = m->getMesh();

        const auto kElastic = m->getMMembrane()->getKElastic();
        const auto eqArea = m->getMMembrane()->getEqArea();

        double area = 0.0;
        for(const auto& t : mesh.getTriangles()) area += t.attr.gTriangle.area;

        const size_t numTriangles = mesh.getTriangles().size();
        for(size_t ti = 0; ti < numTriangles; ++ti) {
            mesh.forEachHalfEdgeInTriangle(ti, [this, &mesh, area, kElastic, eqArea](size_t hei) {
                const auto& dArea = mesh.getHalfEdgeAttribute(hei).gHalfEdge.dTriangleArea;
                _FFType.forces(mesh.getVertexAttribute(mesh.target(hei)).vertex, area, dArea, kElastic, eqArea);
            });
        }
    }
}
