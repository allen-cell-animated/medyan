#ifndef MEDYAN_Structure_SurfaceMesh_Vertex_Hpp
#define MEDYAN_Structure_SurfaceMesh_Vertex_Hpp

#include <memory> // unique_ptr

#include "common.h"
#include "Structure/Bead.h"
#include "Structure/SurfaceMesh/Edge.hpp"
#include "Structure/SurfaceMesh/MVoronoiCell.h"
#include "Structure/SurfaceMesh/Triangle.hpp"

// Forward declarations
class Edge;
class Triangle;
class Membrane;

/******************************************************************************

The vertex class extends bead,
but it is exclusively used in 2D surface meshwork, and contains
information of its neighbors.

******************************************************************************/

class Vertex:
    public Bead, // Inherited from bead to receive full features like coordinate and forces.
                 // But note that in terms of tracking, when the vertex is added to the system,
                 // the base class Bead should also be added to its own collection,
                 // i.e. both the bead and the vertex collection should both have the collection.
    public Database< Vertex, false >
    {
    
    friend class Membrane; // Membrane class can manage id of this vertex

private:
    std::size_t topoIndex_; // Index in the meshwork topology.

    std::unique_ptr<MVoronoiCell> mVertex_; // pointer to Voronoi cell mechanical information

public:
    using vertex_db_type = Database< Vertex, false >;

    ///Main constructor
    Vertex(const Bead::coordinate_type& v, Composite* parent, size_t topoIndex)
        : Bead(mathfunc::vec2Vector(v), parent, 0), topoIndex_(topoIndex)
    {
#ifdef MECHANICS
        // eqArea cannot be obtained at this moment
        mVertex_ = std::make_unique<MVoronoiCell>(getType());
#endif

        usage = Bead::BeadUsage::Membrane;
    }


    void setTopoIndex(size_t index) { topoIndex_ = index; }
    
    // Get mech Voronoi cell
    MVoronoiCell* getMVoronoiCell() { return mVertex_.get(); }

    /// Get all instances of this class from the SubSystem
    static const auto& getVertices() {
        return vertex_db_type::getElements();
    }
    static size_t numVertices() {
        return vertex_db_type::getElements().size();
    }

};

#endif
