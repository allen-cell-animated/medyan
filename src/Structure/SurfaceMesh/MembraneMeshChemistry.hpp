#ifndef MEDYAN_Structure_SurfaceMesh_MembraneMeshChemistry_hpp
#define MEDYAN_Structure_SurfaceMesh_MembraneMeshChemistry_hpp

#include <algorithm>
#include <vector>

#include "Structure/SurfaceMesh/MembraneMeshAttribute.hpp"

namespace medyan {

struct MembraneMeshChemistryInfo {
    // Note:
    //   - indices in this struct must correspond to the actual index as is in
    //     the diffusing species names vector

    struct DiffusionInfo {
        unsigned speciesIndex = 0;

        // Diffusion coeffecient, with dimension L^2 T^-1
        double   diffusionCoeff = 0;
    };
    struct InternalReactionInfo {
        std::vector< unsigned > reactantSpeciesIndices;
        std::vector< unsigned > productSpeciesIndices;

        // Rate constant, with dimension L^(2*(numReactants - 1)) T^-1
        double rateConstant = 0;
    };

    // diffusing species
    std::vector< std::string > diffusingSpeciesNames;

    // diffusion reactions
    std::vector< DiffusionInfo > diffusion;

    // Internal reactions involving species
    std::vector< InternalReactionInfo > internalReactions;
};

inline void setSpeciesAndReactions(
    MembraneMeshAttribute::MeshType&  mesh,
    const MembraneMeshChemistryInfo&  info
) {
    using namespace std;
    using MT = MembraneMeshAttribute::MeshType;

    // Set species
    //---------------------------------
    for(MT::VertexIndex vi {0}; vi < mesh.numVertices(); ++vi) {
        auto& vs = mesh.attribute(vi).vertex->cVertex.species;

        vs.clear();
        for(const auto& name : info.diffusingSpeciesNames) {
            vs.addSpecies(name, 0, max_ulim, RSpeciesType::REG);
        }
    }

    // Set reactions
    //---------------------------------

    const auto indexToSpecies = [&](
        CVertex&                       cv,
        const std::vector< unsigned >& iri
    ) {
        vector< Species* > vs(iri.size());
        for(unsigned i = 0; i < iri.size(); ++i) {
            vs[i] = cv.species.findSpeciesByIndex(i);
        }
        return vs;
    };

    // General reactions
    for(MT::VertexIndex vi {0}; vi < mesh.numVertices(); ++vi) {
        auto& cv = mesh.attribute(vi).vertex->cVertex;

        cv.reactions.clear();
        for(const auto& ri : info.internalReactions) {
            cv.reactions.push_back(make_unique< ReactionDy >(
                indexToSpecies(cv, ri.reactantSpeciesIndices),
                indexToSpecies(cv, ri.productSpeciesIndices),
                ri.rateConstant,
                // volume is not set here
                1.0,
                // rate = rateConstant * area^(1 - numReactants)
                1 - static_cast<int>(ri.reactantSpeciesIndices.size())
            ));
        }
    }

    // Diffusion reactions
    for(MT::HalfEdgeIndex hei {0}; hei < mesh.numHalfEdges(); ++hei) {
        auto& ch     = mesh.attribute(hei).halfEdge->cHalfEdge;
        auto& cvTo   = mesh.attribute(mesh.target(hei)).vertex->cVertex;
        auto& cvFrom = mesh.attribute(mesh.target(mesh.opposite(hei))).vertex->cVertex;

        ch.diffusionReactions.clear();
        for(const auto& di : info.diffusion) {
            ch.diffusionReactions.push_back(make_unique< ReactionDy >(
                indexToSpecies(cvFrom, { di.speciesIndex }),
                indexToSpecies(cvTo,   { di.speciesIndex }),
                di.diffusionCoeff,
                // volume is not set here
                1.0,
                // rate = diffCoeff * area^-1 * shapeFactor
                // where shapeFactor = 0.5 * (cot α + cot β)
                -1
            ));
        }
    }
}


} // namespace medyan

#endif
