/*

Implementing the "easy" bending force field:

TODO: Implement this

*/

#include "MembraneBending.h"

#include "Membrane.h"
#include "Vertex.h"
#include "MVoronoiCell.h"

#include "MembraneBendingVoronoiHelfrich.h"

// Using the Helfrich Hamiltonian of mean curvature in Voronoi cells
double MembraneBending<MembraneBendingVoronoiHelfrich>::computeEnergy(double d) {
    double U = 0;
    double U_i;

    for(auto m: Membrane::getMembranes()) {
        U_i = 0;

        if(d == 0.0) {
            for(Vertex* v: m->getVertexVector()) {
                double kBending = v->getMVoronoiCell()->getBendingModulus();
                double eqCurv = v->getMVoronoiCell()->getEqCurv();

                // The calculation requires that the current area and mean curvature have already been calculated
                double area = v->getMVoronoiCell()->getArea();
                double curv = v->getMVoronoiCell()->getCurv();

                U_i += _FFType.energy(area, curv, kBending, eqCurv); 
            }
        } else {
            for(Vertex *v: m->getVertexVector()) {
                double kBending = v->getMVoronoiCell()->getBendingModulus();
                double eqCurv = v->getMVoronoiCell()->getEqCurv();

                // The calculation requires that the current stretched area and mean curvature have already been calculated
                // As a result, d is just a dummy variable due to its universality
                double areaStretched = v->getMVoronoiCell()->getStretchedArea();
                double curvStretched = v->getMVoronoiCell()->getStretchedCurv();

                U_i += _FFType.energy(areaStretched, curvStretched, kBending, eqCurv, d);
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

void MembraneBending<MembraneBendingVoronoiHelfrich>::computeForces() {
    
    for (auto m: Membrane::getMembranes()) {
    
        for(Vertex* v : m->getVertexVector()){
            
            auto mvc = v->getMVoronoiCell();

            double kBending = mvc->getBendingModulus();
            double eqCurv = mvc->getEqCurv();

            _FFType.forces(v, v->getNeighborVertices(),
                mvc->getArea(), mvc->getDArea(), mvc->getDNeighborArea(),
                mvc->getCurv(), mvc->getDCurv(), mvc->getDNeighborCurv(),
                kBending, eqCurv);
        }
    }
}

void MembraneBending<MembraneBendingVoronoiHelfrich>::computeForcesAux() {
    
    for (auto m: Membrane::getMembranes()) {
    
        for(Vertex* v : m->getVertexVector()){
            
            auto mvc = v->getMVoronoiCell();

            double kBending = mvc->getBendingModulus();
            double eqCurv = mvc->getEqCurv();
           
            _FFType.forcesAux(v, v->getNeighborVertices(),
                mvc->getArea(), mvc->getDArea(), mvc->getDNeighborArea(),
                mvc->getCurv(), mvc->getDCurv(), mvc->getDNeighborCurv(),
                kBending, eqCurv);
        }
    }
}
