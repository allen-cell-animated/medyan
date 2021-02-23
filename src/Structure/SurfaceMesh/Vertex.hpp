#ifndef MEDYAN_Structure_SurfaceMesh_Vertex_hpp
#define MEDYAN_Structure_SurfaceMesh_Vertex_hpp

#include <memory> // unique_ptr
#include <vector>

#include "Chemistry/ReactionDy.hpp"
#include "Chemistry/SpeciesContainer.h"
#include "MathFunctions.h"
#include "Structure/Database.h"
#include "Structure/SurfaceMesh/MVoronoiCell.h"

// CVertex represents a cell around a vertex related to chemistry, and owns
//   - all the diffusing species in the cell
//   - all reactions with only diffusing species in this cell
//
// Note:
//   - diffusion reactions between the cells will be stored in CHalfEdge
struct CVertex {
    using ReactionContainer = std::vector< std::unique_ptr< ReactionDy > >;

    SpeciesPtrContainerVector species;
    ReactionContainer         reactions;
};


class Vertex :
    public Database< Vertex, false >
{
    
private:
    std::size_t topoIndex_; // Index in the meshwork topology.

public:
    using DatabaseType = Database< Vertex, false >;
    using CoordinateType = mathfunc::Vec< 3, floatingpoint >;

    ///Main constructor
    Vertex(const CoordinateType& v, size_t topoIndex)
        : coord(v), force{}, topoIndex_(topoIndex)
    {}


    void setTopoIndex(size_t index) { topoIndex_ = index; }
    
    /// Get all instances of this class from the SubSystem
    static const auto& getVertices() {
        return DatabaseType::getElements();
    }
    static size_t numVertices() {
        return DatabaseType::getElements().size();
    }

    CoordinateType coord;
    CoordinateType force;

    MVertex mVertex; // vertex mechanical information
    CVertex cVertex; // vertex chemical information

};

#endif
