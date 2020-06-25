#ifndef MEDYAN_Structure_SurfaceMesh_MembraneMeshModifier_hpp
#define MEDYAN_Structure_SurfaceMesh_MembraneMeshModifier_hpp

// NOTE
// This file uses both the membrane geometry, as well as the membrane chemistry

#include "Structure/SurfaceMesh/MembraneMeshChemistry.hpp"
#include "Structure/SurfaceMesh/MembraneMeshGeometry.hpp"

namespace medyan {

//-------------------------------------------------------------------------
// Mesh modifiers
//-------------------------------------------------------------------------

// The HalfEdgeMesh already provides basic mesh operations which handles
// the mesh element connections.
// In these operations, some attributes of the mesh need to be handled
// specifically. For example, in material coordinate system, the
// equilibrium area of neighboring triangles after a edge flip need to be
// re-distributed.
// However, these operations cannot be simply achieved by injecting
// behavior (such as passing functions) to the mesh operators of the mesh
// connection class (HalfEdgeMesh), because only with full attribute info
// can one decide what to do/access/modify before/after the mesh
// reconnection.

inline auto insertVertexOnEdge(
    MembraneMeshAttribute::MeshType&             mesh,
    MembraneMeshAttribute::MeshType::EdgeIndex   ei,
    const MembraneMeshAttribute::CoordinateType& newPos
) {
    using MT = MembraneMeshAttribute::MeshType;

    const auto ohei       = mesh.halfEdge(ei);
    const auto ohei_o     = mesh.opposite(ohei);

    const auto opt0       = mesh.polygonType(ohei);
    const bool ist0       = opt0 == MT::PolygonType::triangle;
    const auto opt2       = mesh.polygonType(ohei_o);
    const bool ist2       = opt2 == MT::PolygonType::triangle;

    const auto vi0        = mesh.target(ohei);
    const auto vi2        = mesh.target(ohei_o);

    // Do the vertex insertion
    //---------------------------------
    const auto change = MT::VertexInsertionOnEdge{}(mesh, ei);

    const auto redisEqArea = [](
        double& eqArea1, double a1,
        double& eqArea2, double a2,
        double totalEqArea
    ) {
        // Postcondition: sum of eq areas is equal to the total eq area.
        eqArea1 = totalEqArea * a1 / (a1 + a2);
        eqArea2 = totalEqArea * a2 / (a1 + a2);
    };

    // Set attributes
    //---------------------------------
    mesh.attribute(change.viNew).vertex->coord = newPos;

    // Redistribute properties
    //---------------------------------
    if(ist0) {
        const double a1 = area(mesh, change.tiNew[0]);
        const double a2 = area(mesh, change.tiNew[1]);

        auto& eqArea1 = mesh.attribute(change.tiNew[0]).triangle->mTriangle.eqArea;
        auto& eqArea2 = mesh.attribute(change.tiNew[1]).triangle->mTriangle.eqArea;
        redisEqArea(eqArea1, a1, eqArea2, a2, eqArea1);
    }
    if(ist2) {
        const double a1 = area(mesh, change.tiNew[2]);
        const double a2 = area(mesh, change.tiNew[3]);

        auto& eqArea1 = mesh.attribute(change.tiNew[2]).triangle->mTriangle.eqArea;
        auto& eqArea2 = mesh.attribute(change.tiNew[3]).triangle->mTriangle.eqArea;
        redisEqArea(eqArea1, a1, eqArea2, a2, eqArea1);
    }

    // Relink species and reactions
    // Currently, all species in the new vertex has 0 copy number
    setSpeciesForVertex(mesh, change.viNew, mesh.metaAttribute().chemInfo);
    setInternalReactionsForVertex(mesh, change.viNew, mesh.metaAttribute().chemInfo);
    mesh.forEachHalfEdgeTargetingVertex(change.viNew, [&](MT::HalfEdgeIndex nhei) {
        setDiffusionForHalfEdge(mesh, nhei, mesh.metaAttribute().chemInfo);
        setDiffusionForHalfEdge(mesh, mesh.opposite(nhei), mesh.metaAttribute().chemInfo);
    });

    return change;
}

inline auto collapseEdge(
    MembraneMeshAttribute::MeshType&               mesh,
    MembraneMeshAttribute::MeshType::EdgeIndex     ei,
    std::optional< MembraneMeshAttribute::CoordinateType > newPos
) {
    using namespace std;
    using MT = MembraneMeshAttribute::MeshType;

    const auto hei = mesh.halfEdge(ei);
    const auto hei_o = mesh.opposite(hei);
    const auto vi0 = mesh.target(hei);
    const auto vi1 = mesh.target(hei_o);

    // Record properties
    //---------------------------------
    double oldTotalEqArea = 0.0;
    // Accumulate around two vertices, and subtract triangles on the middle edge
    const auto addOldTotalEqArea = [&](MT::HalfEdgeIndex hei) {
        if(mesh.isInTriangle(hei)) {
            oldTotalEqArea += mesh.attribute(mesh.triangle(hei)).triangle->mTriangle.eqArea;
        }
    };
    const auto subtractOldTotalEqArea = [&](MT::HalfEdgeIndex hei) {
        if(mesh.isInTriangle(hei)) {
            oldTotalEqArea -= mesh.attribute(mesh.triangle(hei)).triangle->mTriangle.eqArea;
        }
    };
    mesh.forEachHalfEdgeTargetingVertex(vi0, addOldTotalEqArea);
    mesh.forEachHalfEdgeTargetingVertex(vi1, addOldTotalEqArea);
    mesh.forEachHalfEdgeInEdge(ei, subtractOldTotalEqArea);

    // Move species vector out of vertices
    auto spe0 = move(mesh.attribute(vi0).vertex->cVertex.species);
    auto spe1 = move(mesh.attribute(vi1).vertex->cVertex.species);
    for(unsigned i = 0; i < spe1.size(); ++i) {
        auto& rs0 = spe0.findSpeciesByIndex(i)->getRSpecies();
        auto& rs1 = spe1.findSpeciesByIndex(i)->getRSpecies();
        rs0.setN(
            rs0.getN() + rs1.getN()
        );
    }
    // spe0 now has the sum of copy numbers.
    // spe1 will be abandoned and destroyed.
    // spe1 cannot be destroyed right now, because there are still reactions
    // pointing to it.

    // Do the edge collapse
    //---------------------------------
    const auto change = MT::EdgeCollapse{}(mesh, ei);

    // Set attributes
    //---------------------------------
    if(newPos.has_value()) mesh.attribute(change.viTo).vertex->coord = *newPos;

    // Redistribute properties
    //---------------------------------
    double newTotalArea = 0.0;
    mesh.forEachHalfEdgeTargetingVertex(change.viTo, [&](MT::HalfEdgeIndex nhei) {
        if(mesh.isInTriangle(nhei)) {
            newTotalArea += area(mesh, mesh.triangle(nhei));
        }
    });
    mesh.forEachHalfEdgeTargetingVertex(change.viTo, [&](MT::HalfEdgeIndex nhei) {
        if(mesh.isInTriangle(nhei)) {
            mesh.attribute(mesh.triangle(nhei)).triangle->mTriangle.eqArea
                = oldTotalEqArea * area(mesh, mesh.triangle(nhei)) / newTotalArea;
        }
    });

    // Relink species and reactions
    mesh.attribute(change.viTo).vertex->cVertex.species = move(spe0);
    setInternalReactionsForVertex(mesh, change.viTo, mesh.metaAttribute().chemInfo);
    mesh.forEachHalfEdgeTargetingVertex(change.viTo, [&](MT::HalfEdgeIndex nhei) {
        setDiffusionForHalfEdge(mesh, nhei, mesh.metaAttribute().chemInfo);
        setDiffusionForHalfEdge(mesh, mesh.opposite(nhei), mesh.metaAttribute().chemInfo);
    });
    // Now spe0 is empty.
    // Now spe1 can be safely destroyed.

    return change;
}

inline void flipEdge(
    MembraneMeshAttribute::MeshType&               mesh,
    MembraneMeshAttribute::MeshType::EdgeIndex     ei
) {
    using MT = MembraneMeshAttribute::MeshType;

    // Precondition:
    //   - ei is not on the border

    // Record old attributes
    //---------------------------------
    double oldTotalEqArea = 0.0;
    mesh.forEachHalfEdgeInEdge(ei, [&](MT::HalfEdgeIndex hei) {
        oldTotalEqArea += mesh.attribute(mesh.triangle(hei)).triangle->mTriangle.eqArea;
    });

    // Do edge flip
    //---------------------------------
    MT::EdgeFlip{}(mesh, ei);

    // Redistribute attributes
    //---------------------------------
    double newTotalArea = 0.0;
    mesh.forEachHalfEdgeInEdge(ei, [&](MT::HalfEdgeIndex hei) {
        newTotalArea += area(mesh, mesh.triangle(hei));
    });
    mesh.forEachHalfEdgeInEdge(ei, [&](MT::HalfEdgeIndex hei) {
        mesh.attribute(mesh.triangle(hei)).triangle->mTriangle.eqArea
            = oldTotalEqArea * area(mesh, mesh.triangle(hei)) / newTotalArea;
    });

    // Relink reactions
    mesh.forEachHalfEdgeInEdge(ei, [&](MT::HalfEdgeIndex hei) {
        setDiffusionForHalfEdge(mesh, hei, mesh.metaAttribute().chemInfo);
    });

}


} // namespace medyan

#endif
