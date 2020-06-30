#ifndef MEDYAN_Structure_SurfaceMesh_HalfEdge_hpp
#define MEDYAN_Structure_SurfaceMesh_HalfEdge_hpp

#include "Chemistry/ReactionDy.hpp"

namespace medyan {

struct CHalfEdge {
    using ReactionContainer = std::vector< std::unique_ptr< ReactionDy > >;

    ReactionContainer diffusionReactions;
};

struct HalfEdge {
    CHalfEdge cHalfEdge;
};

} // namespace medyan

#endif
