#ifndef MEDYAN_Structure_SurfaceMesh_MembraneMeshChemistry_hpp
#define MEDYAN_Structure_SurfaceMesh_MembraneMeshChemistry_hpp

#include <algorithm>
#include <vector>

#include "Structure/SurfaceMesh/MembraneMeshAttribute.hpp"

namespace medyan {

inline void setSpeciesForVertex(
    MembraneMeshAttribute::MeshType&             mesh,
    MembraneMeshAttribute::MeshType::VertexIndex vi,
    const MembraneMeshChemistryInfo&             info
) {
    auto& vs = mesh.attribute(vi).vertex->cVertex.species;

    vs.clear();
    for(const auto& name : info.diffusingSpeciesNames) {
        vs.addSpecies< Species >(name, 0, max_ulim, RSpeciesType::REG);
    }
}

// Helper function: get certain species in CVertex from the species indices
inline auto indicesToSpecies(
    CVertex&                       cv,
    const std::vector< unsigned >& iri
) {
    vector< Species* > vs(iri.size());
    for(unsigned i = 0; i < iri.size(); ++i) {
        vs[i] = cv.species.findSpeciesByIndex(i);
    }
    return vs;
}

// Precondition: Species must have been set in the vertex
inline void setInternalReactionsForVertex(
    MembraneMeshAttribute::MeshType&             mesh,
    MembraneMeshAttribute::MeshType::VertexIndex vi,
    const MembraneMeshChemistryInfo&             info
) {
    auto& cv = mesh.attribute(vi).vertex->cVertex;

    cv.reactions.clear();
    for(const auto& ri : info.internalReactions) {
        cv.reactions.push_back(make_unique< ReactionDy >(
            indicesToSpecies(cv, ri.reactantSpeciesIndices),
            indicesToSpecies(cv, ri.productSpeciesIndices),
            ri.rateConstant,
            // volume is not set here
            1.0,
            // rate = rateConstant * area^(1 - numReactants)
            1 - static_cast<int>(ri.reactantSpeciesIndices.size())
        ));
    }
}

// Precondition: Species must have been set in the vertices
inline void setDiffusionForHalfEdge(
    MembraneMeshAttribute::MeshType&               mesh,
    MembraneMeshAttribute::MeshType::HalfEdgeIndex hei,
    const MembraneMeshChemistryInfo&               info
) {
    auto& ch     = mesh.attribute(hei).halfEdge->cHalfEdge;
    auto& cvTo   = mesh.attribute(mesh.target(hei)).vertex->cVertex;
    auto& cvFrom = mesh.attribute(mesh.target(mesh.opposite(hei))).vertex->cVertex;

    ch.diffusionReactions.clear();
    for(const auto& di : info.diffusion) {
        ch.diffusionReactions.push_back(make_unique< ReactionDy >(
            indicesToSpecies(cvFrom, { di.speciesIndex }),
            indicesToSpecies(cvTo,   { di.speciesIndex }),
            di.diffusionCoeff,
            // volume is not set here
            1.0,
            // rate = diffCoeff * area^-1 * shapeFactor
            // where shapeFactor = 0.5 * (cot α + cot β)
            -1
        ));
    }
}

inline void setSpeciesAndReactions(
    MembraneMeshAttribute::MeshType&  mesh,
    const MembraneMeshChemistryInfo&  info
) {
    using namespace std;
    using MT = MembraneMeshAttribute::MeshType;

    // Set species
    //---------------------------------
    for(MT::VertexIndex vi {0}; vi < mesh.numVertices(); ++vi) {
        setSpeciesForVertex(mesh, vi, info);
    }

    // Set reactions
    //---------------------------------

    // General reactions
    for(MT::VertexIndex vi {0}; vi < mesh.numVertices(); ++vi) {
        setInternalReactionsForVertex(mesh, vi, info);
    }

    // Diffusion reactions
    for(MT::HalfEdgeIndex hei {0}; hei < mesh.numHalfEdges(); ++hei) {
        setDiffusionForHalfEdge(mesh, hei, info);
    }
}


} // namespace medyan

#endif
