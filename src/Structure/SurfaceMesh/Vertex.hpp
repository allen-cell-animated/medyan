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
    
private:
    std::size_t topoIndex_; // Index in the meshwork topology.

    std::unique_ptr<MVoronoiCell> mVertex_; // pointer to Voronoi cell mechanical information

public:
    using DatabaseType = Database< Vertex, false >;
    using CoordinateType = mathfunc::Vec< 3, floatingpoint >;

    ///Main constructor
    Vertex(const CoordinateType& v, size_t topoIndex)
        : coord(v), force{}, topoIndex_(topoIndex)
    {
        // eqArea cannot be obtained at this moment
        mVertex_ = std::make_unique<MVoronoiCell>(getType());

    }


    void setTopoIndex(size_t index) { topoIndex_ = index; }
    
    // Get mech Voronoi cell
    MVoronoiCell* getMVoronoiCell() { return mVertex_.get(); }

    /// Get all instances of this class from the SubSystem
    static const auto& getVertices() {
        return DatabaseType::getElements();
    }
    static size_t numVertices() {
        return DatabaseType::getElements().size();
    }

    CoordinateType coord;
    CoordinateType force;

};

#endif
