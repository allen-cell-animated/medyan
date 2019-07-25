
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

#include "Mechanics/ForceField/Volume/TriangleCylinderExclVolume.hpp"

#include <algorithm> // max

#include "MathFunctions.h"
using namespace mathfunc;

#include "Mechanics/ForceField/Volume/TriangleCylinderBeadExclVolRepulsion.hpp"

#include "Structure/Bead.h"
#include "Structure/Cylinder.h"
#include "Structure/SurfaceMesh/Membrane.hpp"
#include "Structure/SurfaceMesh/Triangle.hpp"
#include "Structure/SurfaceMesh/Vertex.h"
#include "Util/Math/RayTriangleIntersect.hpp"

template <class TriangleCylinderExclVolumeInteractionType>
floatingpoint TriangleCylinderExclVolume<TriangleCylinderExclVolumeInteractionType>::computeEnergy(const floatingpoint* coord, bool stretched) {
    
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
                    static_cast<Vec3>(makeVec<3>(coord + 3 * v0->Bead::getStableIndex())),
                    static_cast<Vec3>(makeVec<3>(coord + 3 * v1->Bead::getStableIndex())),
                    static_cast<Vec3>(makeVec<3>(coord + 3 * v2->Bead::getStableIndex())),
                    static_cast<Vec3>(makeVec<3>(coord + 3 * b->getStableIndex())),
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
void TriangleCylinderExclVolume<TriangleCylinderExclVolumeInteractionType>::computeForces(const floatingpoint* coord, floatingpoint* force) {

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
                    force + 3 * v0->Bead::getStableIndex(),
                    force + 3 * v1->Bead::getStableIndex(),
                    force + 3 * v2->Bead::getStableIndex(),
                    force + 3 * b->getStableIndex(),
                    static_cast<Vec3>(makeVec<3>(coord + 3 * v0->Bead::getStableIndex())),
                    static_cast<Vec3>(makeVec<3>(coord + 3 * v1->Bead::getStableIndex())),
                    static_cast<Vec3>(makeVec<3>(coord + 3 * v2->Bead::getStableIndex())),
                    static_cast<Vec3>(makeVec<3>(coord + 3 * b->getStableIndex())),
                    area, dArea0, dArea1, dArea2, kExVol);
            }
        }
    }
}

namespace {

template< typename InteractionType, typename VecType >
void exclVolLoadForce(
    const InteractionType& interaction, double area, double kExVol,
    const Bead& bo, Bead& bd, typename TriangleCylinderExclVolume< InteractionType >::LoadForceEnd end,
    const VecType& v0, const VecType& v1, const VecType& v2
) {
    using LoadForceEnd = typename TriangleCylinderExclVolume< InteractionType >::LoadForceEnd;

    auto& loadForces = (end == LoadForceEnd::Plus ? bd.loadForcesP : bd.loadForcesM);
    auto& lfi        = (end == LoadForceEnd::Plus ? bd.lfip        : bd.lfim       );

    // This is the direction of polymerization
    const auto dir = normalizedVector(bd.coordinate() - bo.coordinate());

    // Test intersection of the ray with the triangle
    const auto intersectRes = ray_tracing::MollerTrumboreIntersect<>()(
        bd.coordinate(), dir,
        v0, v1, v2
    );

    //array of coordinate values to update
    const auto monSize = SysParams::Geometry().monomerSize   [bd.getType()];
    auto       cylSize = SysParams::Geometry().cylinderNumMon[bd.getType()];

    // Set load force to inf right before intersection.
    if(intersectRes.intersect && intersectRes.t > 0) {
        const int maxCylSize = static_cast<int>(intersectRes.t / monSize);
        if(maxCylSize < cylSize) {
            loadForces[maxCylSize] = std::numeric_limits< floatingpoint >::infinity();
            cylSize = maxCylSize;
        }
    }

    for (int i = 0; i < cylSize; i++) {

        const auto newCoord = bd.coordinate() + (i * monSize) * dir;

        const auto loadForce = interaction.loadForces(v0, v1, v2, newCoord, area, kExVol);
        const auto effLoadForce = std::max< double >(-dot(dir, loadForce), 0.0);

        loadForces[i] += effLoadForce;
    }
    //reset lfi
    lfi = 0;

} // void exclVolLoadForce(...)

} // namespace (anonymous)

template< typename InteractionType >
void TriangleCylinderExclVolume< InteractionType >::computeLoadForce(const Bead& bo, Bead& bd, LoadForceEnd end, const Triangle& t) const {
    const auto& mesh = t.getParent()->getMesh();
    const size_t ti = t.getTopoIndex();
    const size_t hei0 = mesh.getTriangles()[ti].halfEdgeIndex;
    const size_t hei1 = mesh.next(hei0);
    const size_t hei2 = mesh.next(hei1);
    const Vec3 v0 (mesh.getVertexAttribute(mesh.target(hei0)).getCoordinate());
    const Vec3 v1 (mesh.getVertexAttribute(mesh.target(hei1)).getCoordinate());
    const Vec3 v2 (mesh.getVertexAttribute(mesh.target(hei2)).getCoordinate());

    const auto area = mesh.getTriangleAttribute(ti).gTriangle.area;
    double kExVol = t.getMTriangle()->getExVolConst();

    exclVolLoadForce(
        _FFType, area, kExVol,
        bo, bd, end,
        v0, v1, v2
    );
}

template < typename InteractionType >
void TriangleCylinderExclVolume< InteractionType >::computeLoadForces() {

    for(auto t: Triangle::getTriangles()) {

        const auto& mesh = t->getParent()->getMesh();
        const size_t ti = t->getTopoIndex();
        const size_t hei0 = mesh.getTriangles()[ti].halfEdgeIndex;
        const size_t hei1 = mesh.next(hei0);
        const size_t hei2 = mesh.next(hei1);
        const Vec3 v0 (mesh.getVertexAttribute(mesh.target(hei0)).getCoordinate());
        const Vec3 v1 (mesh.getVertexAttribute(mesh.target(hei1)).getCoordinate());
        const Vec3 v2 (mesh.getVertexAttribute(mesh.target(hei2)).getCoordinate());

        const auto area = mesh.getTriangleAttribute(ti).gTriangle.area;
        double kExVol = t->getMTriangle()->getExVolConst();
        
        for(auto &c: _neighborList->getNeighbors(t)) {

            // potential acts on second cylinder bead if it is a plus  end
            // potential acts on first  cylinder bead if it is a minus end
            if(c->isPlusEnd()) {
                exclVolLoadForce(
                    _FFType, area, kExVol,
                    *c->getFirstBead(), *c->getSecondBead(), LoadForceEnd::Plus,
                    v0, v1, v2
                );
            }
            
            if(c->isMinusEnd()) {
                exclVolLoadForce(
                    _FFType, area, kExVol,
                    *c->getSecondBead(), *c->getFirstBead(), LoadForceEnd::Minus,
                    v0, v1, v2
                );
            }
            
        }
    }

} // void ...::computeLoadForces()

// Template instantiation
template class TriangleCylinderExclVolume< TriangleCylinderBeadExclVolRepulsion >;
