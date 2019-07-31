
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

#include "Mechanics/ForceField/Volume/TriangleBeadExclVolume.hpp"

#include <algorithm> // max

#include "MathFunctions.h"
using namespace mathfunc;

#include "Mechanics/ForceField/Volume/TriangleBeadExclVolRepulsion.hpp"

#include "Structure/Bead.h"
#include "Structure/Cylinder.h"
#include "Structure/Filament.h"
#include "Structure/SurfaceMesh/Membrane.hpp"
#include "Structure/SurfaceMesh/Triangle.hpp"
#include "Structure/SurfaceMesh/Vertex.hpp"
#include "Util/Math/RayTriangleIntersect.hpp"

template< typename InteractionType >
floatingpoint TriangleBeadExclVolume< InteractionType >::computeEnergy(const floatingpoint* coord, bool stretched) {
    
    double U = 0;
    double U_i;

    for(auto t: Triangle::getTriangles()) {

        auto& mesh = t->getParent()->getMesh();
        Membrane::MembraneMeshAttributeType::cacheIndices(mesh);

        const size_t ti = t->getTopoIndex();
        const auto& ta = mesh.getTriangleAttribute(ti);
        const size_t vi0 = ta.cachedCoordIndex[0];
        const size_t vi1 = ta.cachedCoordIndex[1];
        const size_t vi2 = ta.cachedCoordIndex[2];

        const auto area = stretched ? ta.gTriangleS.area : ta.gTriangle.area;
        const double kExVol = t->getMTriangle()->getExVolConst();
        
        if(_neighborList->hasNeighbor(t)) for(auto b : _neighborList->getNeighbors(t)) {
                
            U_i = _FFType.energy(
                static_cast<Vec3d>(makeVec<3>(coord + 3 * vi0)),
                static_cast<Vec3d>(makeVec<3>(coord + 3 * vi1)),
                static_cast<Vec3d>(makeVec<3>(coord + 3 * vi2)),
                static_cast<Vec3d>(makeVec<3>(coord + 3 * b->getStableIndex())),
                area, kExVol);
            
            if(!std::isfinite(U_i) || U_i < -1.0) {
                
                //set culprits and exit
                triangleCulprit_ = t;
                beadCulprit_     = b;
                
                return -1;
            }
            else
                U += U_i;
        }
    }
    
    return U;
}

template< typename InteractionType >
void TriangleBeadExclVolume< InteractionType >::computeForces(const floatingpoint* coord, floatingpoint* force) {

    for(auto t: Triangle::getTriangles()) {

        auto& mesh = t->getParent()->getMesh();
        Membrane::MembraneMeshAttributeType::cacheIndices(mesh);

        const size_t ti = t->getTopoIndex();
        const auto& ta = mesh.getTriangleAttribute(ti);
        const size_t hei0 = ta.cachedHalfEdgeIndex[0];
        const size_t hei1 = ta.cachedHalfEdgeIndex[1];
        const size_t hei2 = ta.cachedHalfEdgeIndex[2];
        const size_t vi0 = ta.cachedCoordIndex[0];
        const size_t vi1 = ta.cachedCoordIndex[1];
        const size_t vi2 = ta.cachedCoordIndex[2];

        const auto area = ta.gTriangle.area;
        const auto& dArea0 = mesh.getHalfEdgeAttribute(hei0).gHalfEdge.dTriangleArea;
        const auto& dArea1 = mesh.getHalfEdgeAttribute(hei1).gHalfEdge.dTriangleArea;
        const auto& dArea2 = mesh.getHalfEdgeAttribute(hei2).gHalfEdge.dTriangleArea;
        const double kExVol = t->getMTriangle()->getExVolConst();

        if(_neighborList->hasNeighbor(t)) for(auto b : _neighborList->getNeighbors(t)) {

            _FFType.forces(
                force + 3 * vi0,
                force + 3 * vi1,
                force + 3 * vi2,
                force + 3 * b->getStableIndex(),
                static_cast<Vec3>(makeVec<3>(coord + 3 * vi0)),
                static_cast<Vec3>(makeVec<3>(coord + 3 * vi1)),
                static_cast<Vec3>(makeVec<3>(coord + 3 * vi2)),
                static_cast<Vec3>(makeVec<3>(coord + 3 * b->getStableIndex())),
                area, dArea0, dArea1, dArea2, kExVol);
        }
    }
}

namespace {

template< typename InteractionType, typename VecType >
void exclVolLoadForce(
    const InteractionType& interaction, double area, double kExVol,
    const Bead& bo, Bead& bd, typename TriangleBeadExclVolume< InteractionType >::LoadForceEnd end,
    const VecType& v0, const VecType& v1, const VecType& v2
) {
    using LoadForceEnd = typename TriangleBeadExclVolume< InteractionType >::LoadForceEnd;

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
void TriangleBeadExclVolume< InteractionType >::computeLoadForces() const {

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
        
        if(_neighborList->hasNeighbor(t)) for(auto b : _neighborList->getNeighbors(t)) {

            // potential acts on second cylinder bead if it is a plus  end
            // potential acts on first  cylinder bead if it is a minus end
            Cylinder& cPlus = *static_cast< Filament* >(b->getParent())->getPlusEndCylinder();
            if(b == cPlus.getSecondBead()) {
                exclVolLoadForce(
                    _FFType, area, kExVol,
                    *cPlus.getFirstBead(), *b, LoadForceEnd::Plus,
                    v0, v1, v2
                );
            }

            Cylinder& cMinus = *static_cast< Filament* >(b->getParent())->getMinusEndCylinder();
            if(b == cMinus.getFirstBead()) {
                exclVolLoadForce(
                    _FFType, area, kExVol,
                    *cMinus.getSecondBead(), *b, LoadForceEnd::Minus,
                    v0, v1, v2
                );
            }
            
        }
    }

} // void ...::computeLoadForces() const

template< typename InteractionType >
void TriangleBeadExclVolume< InteractionType >::computeLoadForce(Cylinder* c, LoadForceEnd end) const {

    Bead* tip   = (end == LoadForceEnd::Plus ? c->getSecondBead() : c->getFirstBead());
    Bead* other = (end == LoadForceEnd::Plus ? c->getFirstBead() : c->getSecondBead());

    if(_neighborList->hasNeighbor(tip)) for(auto t : _neighborList->getNeighbors(tip)) {

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

        exclVolLoadForce(
            _FFType, area, kExVol,
            *other, *tip, end,
            v0, v1, v2
        );
    }

} // void ...::computeLoadForce(...) const

// Template instantiation
template class TriangleBeadExclVolume< TriangleBeadExclVolRepulsion >;
