#ifndef MEDYAN_Structure_SurfaceMesh_Vertex_Hpp
#define MEDYAN_Structure_SurfaceMesh_Vertex_Hpp

#include <memory> // unique_ptr

#include "common.h"

#include "Bead.h"

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
    size_t _topoIndex; // Index in the meshwork topology.

    unique_ptr<MVoronoiCell> _mVoronoiCell; // pointer to Voronoi cell mechanical information

    size_t _membraneVertexIdx; // The index of this in the _vertexVector of the parent membrane

public:
    using vertex_db_type = Database< Vertex, false >;

    ///Main constructor
    Vertex(const Bead::coordinate_type& v, Composite* parent, size_t topoIndex);

    void setTopoIndex(size_t index) { _topoIndex = index; }
    
    // Get mech Voronoi cell
    MVoronoiCell* getMVoronoiCell() { return _mVoronoiCell.get(); }

    /// Get all instances of this class from the SubSystem
    static const auto& getVertices() {
        return vertex_db_type::getElements();
    }
    static size_t numVertices() {
        return vertex_db_type::getElements().size();
    }

    size_t getMembraneVertexIdx()const { return _membraneVertexIdx; }

    //@{
    /// SubSystem management, inherited from Bead: Trackable
    // Incremental to the functions in Bead class.
    virtual void addToSubSystem()override {}
    virtual void removeFromSubSystem()override {}
    //@}

};

#endif
