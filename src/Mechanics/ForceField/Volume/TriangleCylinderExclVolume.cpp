
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

#include "Structure/SurfaceMesh/Triangle.h"
#include "Structure/SurfaceMesh/Vertex.h"
#include "Cylinder.h"
#include "Bead.h"
#include "Structure/SurfaceMesh/Membrane.hpp"
#include "util/math/RayTriangleIntersect.hpp"

template <class TriangleCylinderExclVolumeInteractionType>
double TriangleCylinderExclVolume<TriangleCylinderExclVolumeInteractionType>::computeEnergy(const double* coord, bool stretched) {
    
    double U = 0;
    double U_i;

    for(auto t: Triangle::getTriangles()) {

        const auto& mesh = t->getParent()->getMesh();
        const size_t ti = t->getTopoIndex();
        const size_t hei0 = mesh.getTriangles()[ti].halfEdgeIndex;
        const size_t hei1 = mesh.next(hei0);
        const size_t hei2 = mesh.next(hei1);
        Vertex* const v0 = mesh.getVertexAttribute(mesh.target(hei0)).vertex;
        Vertex* const v1 = mesh.getVertexAttribute(mesh.target(hei1)).vertex;
        Vertex* const v2 = mesh.getVertexAttribute(mesh.target(hei2)).vertex;

        const auto area = stretched ?
            mesh.getTriangleAttribute(ti).gTriangleS.area :
            mesh.getTriangleAttribute(ti).gTriangle.area;
        double kExVol = t->getMTriangle()->getExVolConst();
        
        for(auto &c : _neighborList->getNeighbors(t)) {

            // Use only 1st bead unless the cylinder is at plus end
            size_t numBeads = (c->isPlusEnd()? 2: 1);

            for(size_t idx = 0; idx < numBeads; ++idx) {
                Bead* b = (idx? c->getSecondBead(): c->getFirstBead());
                
                U_i = _FFType.energy(
                    makeVec<3>(coord + 3 * v0->Bead::getIndex()),
                    makeVec<3>(coord + 3 * v1->Bead::getIndex()),
                    makeVec<3>(coord + 3 * v2->Bead::getIndex()),
                    makeVec<3>(coord + 3 * b->getIndex()),
                    area, kExVol);
                
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
void TriangleCylinderExclVolume<TriangleCylinderExclVolumeInteractionType>::computeForces(const double* coord, double* force) {

    for(auto t: Triangle::getTriangles()) {

        const auto& mesh = t->getParent()->getMesh();
        const size_t ti = t->getTopoIndex();
        const size_t hei0 = mesh.getTriangles()[ti].halfEdgeIndex;
        const size_t hei1 = mesh.next(hei0);
        const size_t hei2 = mesh.next(hei1);
        Vertex* const v0 = mesh.getVertexAttribute(mesh.target(hei0)).vertex;
        Vertex* const v1 = mesh.getVertexAttribute(mesh.target(hei1)).vertex;
        Vertex* const v2 = mesh.getVertexAttribute(mesh.target(hei2)).vertex;

        const auto area = mesh.getTriangleAttribute(ti).gTriangle.area;
        const auto& dArea0 = mesh.getHalfEdgeAttribute(hei0).gHalfEdge.dTriangleArea;
        const auto& dArea1 = mesh.getHalfEdgeAttribute(hei1).gHalfEdge.dTriangleArea;
        const auto& dArea2 = mesh.getHalfEdgeAttribute(hei2).gHalfEdge.dTriangleArea;
        double kExVol = t->getMTriangle()->getExVolConst();

        for(auto &c: _neighborList->getNeighbors(t)) {

            // Use only 1st bead unless the cylinder is at plus end
            size_t numBeads = (c->isPlusEnd()? 2: 1);

            for(size_t idx = 0; idx < numBeads; ++idx) {
                Bead* b = (idx? c->getSecondBead(): c->getFirstBead());
            
                _FFType.forces(
                    force + 3 * v0->Bead::getIndex(),
                    force + 3 * v1->Bead::getIndex(),
                    force + 3 * v2->Bead::getIndex(),
                    force + 3 * b->getIndex(),
                    makeVec<3>(coord + 3 * v0->Bead::getIndex()),
                    makeVec<3>(coord + 3 * v1->Bead::getIndex()),
                    makeVec<3>(coord + 3 * v2->Bead::getIndex()),
                    makeVec<3>(coord + 3 * b->getIndex()),
                    area, dArea0, dArea1, dArea2, kExVol);
            }
        }
    }
}

template <class TriangleCylinderExclVolumeInteractionType>
void TriangleCylinderExclVolume<TriangleCylinderExclVolumeInteractionType>::computeLoadForces() {

    for(auto t: Triangle::getTriangles()) {

        const auto& mesh = t->getParent()->getMesh();
        const size_t ti = t->getTopoIndex();
        const size_t hei0 = mesh.getTriangles()[ti].halfEdgeIndex;
        const size_t hei1 = mesh.next(hei0);
        const size_t hei2 = mesh.next(hei1);
        const Vec3 v0 = mesh.getVertexAttribute(mesh.target(hei0)).getCoordinate();
        const Vec3 v1 = mesh.getVertexAttribute(mesh.target(hei1)).getCoordinate();
        const Vec3 v2 = mesh.getVertexAttribute(mesh.target(hei2)).getCoordinate();

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
                const auto normal = normalizedVector(bd->coordinate() - bo->coordinate());

                // Test intersection of the ray with the triangle
                const auto intersectRes = ray_tracing::MollerTrumboreIntersect<>()(
                    bd->coordinate(), normal,
                    v0, v1, v2
                );

                //array of coordinate values to update
                const auto monSize = SysParams::Geometry().monomerSize[bd->getType()];
                auto cylSize = SysParams::Geometry().cylinderNumMon[bd->getType()];

                // Set load force to inf right before intersection.
                if(intersectRes.intersect && intersectRes.t > 0) {
                    const int maxCylSize = static_cast<int>(intersectRes.t / monSize);
                    if(maxCylSize < cylSize) {
                        bd->loadForcesP[maxCylSize] = std::numeric_limits<double>::infinity();
                        cylSize = maxCylSize;
                    }
                }

                bd->lfip = 0;
                for (int i = 0; i < cylSize; i++) {
                    
                    Vec3 newCoord {
                        bd->coordinate()[0] + i * normal[0] * monSize,
                        bd->coordinate()[1] + i * normal[1] * monSize,
                        bd->coordinate()[2] + i * normal[2] * monSize
                    };
                    
                    auto loadForce = _FFType.loadForces(v0, v1, v2, newCoord, area, kExVol); // FIXME change it
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
                const auto normal = normalizedVector(bd->coordinate() - bo->coordinate());

                // Test intersection of the ray with the triangle
                const auto intersectRes = ray_tracing::MollerTrumboreIntersect<>()(
                    bd->coordinate(), normal,
                    v0, v1, v2
                );

                //array of coordinate values to update
                const auto monSize = SysParams::Geometry().monomerSize[bd->getType()];
                auto cylSize = SysParams::Geometry().cylinderNumMon[bd->getType()];

                // Set load force to inf right before intersection.
                if(intersectRes.intersect && intersectRes.t > 0) {
                    const int maxCylSize = static_cast<int>(intersectRes.t / monSize);
                    if(maxCylSize < cylSize) {
                        bd->loadForcesP[maxCylSize] = std::numeric_limits<double>::infinity();
                        cylSize = maxCylSize;
                    }
                }

                bd->lfim = 0;
                for (int i = 0; i < cylSize; i++) {
                    
                    Vec3 newCoord {
                        bd->coordinate()[0] + i * normal[0] * monSize,
                        bd->coordinate()[1] + i * normal[1] * monSize,
                        bd->coordinate()[2] + i * normal[2] * monSize
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
template double TriangleCylinderExclVolume<TriangleCylinderBeadExclVolRepulsion>::computeEnergy(const double* coord, bool stretched);
template void TriangleCylinderExclVolume<TriangleCylinderBeadExclVolRepulsion>::computeForces(const double* coord, double* force);
template void TriangleCylinderExclVolume<TriangleCylinderBeadExclVolRepulsion>::computeLoadForces();
