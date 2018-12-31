
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v3.1
//
//  Copyright (2017-2018)  Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the MEDYAN web page for more information:
//  http://www.medyan.org
//------------------------------------------------------------------

#include "TriangleCylinderExclVolume.h"

#include "MathFunctions.h"
using namespace mathfunc;

#include "TriangleCylinderBeadExclVolRepulsion.h"

#include "Triangle.h"
#include "Vertex.h"
#include "Cylinder.h"
#include "Bead.h"
#include "Structure/SurfaceMesh/Membrane.hpp"

template <class TriangleCylinderExclVolumeInteractionType>
double TriangleCylinderExclVolume<TriangleCylinderExclVolumeInteractionType>::computeEnergy(bool stretched) {
    
    double U = 0;
    double U_i;

    for(auto t: Triangle::getTriangles()) {

        const auto& mesh = static_cast<Membrane*>(t->getParent())->getMesh();
        const size_t ti = t->getTopoIndex();
        const size_t hei0 = mesh.getTriangles()[ti].halfEdgeIndex;
        const size_t hei1 = mesh.next(hei0);
        const size_t hei2 = mesh.next(hei1);
        const Vertex* v0 = mesh.getVertexAttribute(mesh.target(hei0)).vertex;
        const Vertex* v1 = mesh.getVertexAttribute(mesh.target(hei1)).vertex;
        const Vertex* v2 = mesh.getVertexAttribute(mesh.target(hei2)).vertex;

        const auto area = mesh.getTriangleAttribute(ti).gTriangle.area;
        double kExVol = t->getMTriangle()->getExVolConst();
        
        for(auto &c : _neighborList->getNeighbors(t)) {

            // Use only 1st bead unless the cylinder is at plus end
            size_t numBeads = (c->isPlusEnd()? 2: 1);

            for(size_t idx = 0; idx < numBeads; ++idx) {
                Bead* b = (idx? c->getSecondBead(): c->getFirstBead());
                
                U_i = _FFType.energy(v0, v1, v2, b, area, kExVol, stretched);
                
                if(fabs(U_i) == numeric_limits<double>::infinity()
                || U_i != U_i || U_i < -1.0) {
                    
                    //set culprits and exit
                    _triangleCulprit = t;
                    _cylinderCulprit = c;
                    
                    return -1;
                }
                else
                    U += U_i;
            }
        }
    }
    
    return U;
}

template <class TriangleCylinderExclVolumeInteractionType>
void TriangleCylinderExclVolume<TriangleCylinderExclVolumeInteractionType>::computeForces() {

    for(auto t: Triangle::getTriangles()) {

        const auto& mesh = static_cast<Membrane*>(t->getParent())->getMesh();
        const size_t ti = t->getTopoIndex();
        const size_t hei0 = mesh.getTriangles()[ti].halfEdgeIndex;
        const size_t hei1 = mesh.next(hei0);
        const size_t hei2 = mesh.next(hei1);
        Vertex* const v0 = mesh.getVertexAttribute(mesh.target(hei0)).vertex;
        Vertex* const v1 = mesh.getVertexAttribute(mesh.target(hei1)).vertex;
        Vertex* const v2 = mesh.getVertexAttribute(mesh.target(hei2)).vertex;

        const auto area = mesh.getTriangleAttribute(ti).gTriangle.area;
        const auto& dArea0 = mesh.getHalfEdgeAttribute(hei0).gHalfEdge.dArea;
        const auto& dArea1 = mesh.getHalfEdgeAttribute(hei1).gHalfEdge.dArea;
        const auto& dArea2 = mesh.getHalfEdgeAttribute(hei2).gHalfEdge.dArea;
        double kExVol = t->getMTriangle()->getExVolConst();

        for(auto &c: _neighborList->getNeighbors(t)) {

            // Use only 1st bead unless the cylinder is at plus end
            size_t numBeads = (c->isPlusEnd()? 2: 1);

            for(size_t idx = 0; idx < numBeads; ++idx) {
                Bead* b = (idx? c->getSecondBead(): c->getFirstBead());
            
                _FFType.forces(v0, v1, v2, b, area, dArea0, dArea1, dArea2, kExVol);
            }
        }
    }
}


template <class TriangleCylinderExclVolumeInteractionType>
void TriangleCylinderExclVolume<TriangleCylinderExclVolumeInteractionType>::computeForcesAux() {

    for(auto t: Triangle::getTriangles()) {

        const auto& mesh = static_cast<Membrane*>(t->getParent())->getMesh();
        const size_t ti = t->getTopoIndex();
        const size_t hei0 = mesh.getTriangles()[ti].halfEdgeIndex;
        const size_t hei1 = mesh.next(hei0);
        const size_t hei2 = mesh.next(hei1);
        Vertex* const v0 = mesh.getVertexAttribute(mesh.target(hei0)).vertex;
        Vertex* const v1 = mesh.getVertexAttribute(mesh.target(hei1)).vertex;
        Vertex* const v2 = mesh.getVertexAttribute(mesh.target(hei2)).vertex;

        const auto area = mesh.getTriangleAttribute(ti).gTriangle.area;
        const auto& dArea0 = mesh.getHalfEdgeAttribute(hei0).gHalfEdge.dArea;
        const auto& dArea1 = mesh.getHalfEdgeAttribute(hei1).gHalfEdge.dArea;
        const auto& dArea2 = mesh.getHalfEdgeAttribute(hei2).gHalfEdge.dArea;
        double kExVol = t->getMTriangle()->getExVolConst();
        
        for(auto &c: _neighborList->getNeighbors(t)) {

            // Use only 1st bead unless the cylinder is at plus end
            size_t numBeads = (c->isPlusEnd()? 2: 1);

            for(size_t idx = 0; idx < numBeads; ++idx) {
                Bead* b = (idx? c->getSecondBead(): c->getFirstBead());
            
                _FFType.forcesAux(v0, v1, v2, b, area, dArea0, dArea1, dArea2, kExVol);
            }
        }
    }
}

template <class TriangleCylinderExclVolumeInteractionType>
void TriangleCylinderExclVolume<TriangleCylinderExclVolumeInteractionType>::computeLoadForces() {

    for(auto t: Triangle::getTriangles()) {

        const auto& mesh = static_cast<Membrane*>(t->getParent())->getMesh();
        const size_t ti = t->getTopoIndex();
        const size_t hei0 = mesh.getTriangles()[ti].halfEdgeIndex;
        const size_t hei1 = mesh.next(hei0);
        const size_t hei2 = mesh.next(hei1);
        const Vertex* v0 = mesh.getVertexAttribute(mesh.target(hei0)).vertex;
        const Vertex* v1 = mesh.getVertexAttribute(mesh.target(hei1)).vertex;
        const Vertex* v2 = mesh.getVertexAttribute(mesh.target(hei2)).vertex;

        const auto area = mesh.getTriangleAttribute(ti).gTriangle.area;
        double kExVol = t->getMTriangle()->getExVolConst();
        
        for(auto &c: _neighborList->getNeighbors(t)) {

            // potential acts on second cylinder bead if it is a plus  end
            // potential acts on first  cylinder bead if it is a minus end
            Bead* bd;
            Bead* bo;
            if(c->isPlusEnd()) {
                
                bd = c->getSecondBead();
                bo = c->getFirstBead();
                
                ///this normal is in the direction of polymerization
                auto normal = vector2Vec<3, double>(twoPointDirection(bo->coordinate, bd->coordinate));
                
                //array of coordinate values to update
                auto monSize = SysParams::Geometry().monomerSize[bd->getType()];
                auto cylSize = SysParams::Geometry().cylinderNumMon[bd->getType()];
                
                bd->lfip = 0;
                for (int i = 0; i < cylSize; i++) {
                    
                    Vec3 newCoord {
                        bd->coordinate[0] + i * normal[0] * monSize,
                        bd->coordinate[1] + i * normal[1] * monSize,
                        bd->coordinate[2] + i * normal[2] * monSize
                    };
                    
                    auto loadForce = _FFType.loadForces(v0, v1, v2, newCoord, area, kExVol);
                    double effLoadForce = -dot(normal, loadForce);
                    if(effLoadForce < 0.0) effLoadForce = 0.0;

                    bd->loadForcesP[bd->lfip++] += effLoadForce;
                }
                //reset lfi
                bd->lfip = 0;
            }
            
            if(c->isMinusEnd()) {
                
                bd = c->getFirstBead();
                bo = c->getSecondBead();
                
                ///this normal is in the direction of polymerization
                auto normal = vector2Vec<3, double>(twoPointDirection(bo->coordinate, bd->coordinate));
                
                //array of coordinate values to update
                auto monSize = SysParams::Geometry().monomerSize[bd->getType()];
                auto cylSize = SysParams::Geometry().cylinderNumMon[bd->getType()];
                
                
                bd->lfim = 0;
                for (int i = 0; i < cylSize; i++) {
                    
                    Vec3 newCoord {
                        bd->coordinate[0] + i * normal[0] * monSize,
                        bd->coordinate[1] + i * normal[1] * monSize,
                        bd->coordinate[2] + i * normal[2] * monSize
                    };
                    
                    auto loadForce = _FFType.loadForces(v0, v1, v2, newCoord, area, kExVol);
                    double effLoadForce = -dot(normal, loadForce);
                    if(effLoadForce < 0.0) effLoadForce = 0.0;

                    bd->loadForcesM[bd->lfim++] += effLoadForce;
                }
                //reset lfi
                bd->lfim = 0;
            }
            
        }
    }

}

///Template specializations
template double TriangleCylinderExclVolume<TriangleCylinderBeadExclVolRepulsion>::computeEnergy(bool stretched);
template void TriangleCylinderExclVolume<TriangleCylinderBeadExclVolRepulsion>::computeForces();
template void TriangleCylinderExclVolume<TriangleCylinderBeadExclVolRepulsion>::computeForcesAux();
template void TriangleCylinderExclVolume<TriangleCylinderBeadExclVolRepulsion>::computeLoadForces();
