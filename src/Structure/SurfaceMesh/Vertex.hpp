#ifndef MEDYAN_Structure_SurfaceMesh_Vertex_hpp
#define MEDYAN_Structure_SurfaceMesh_Vertex_hpp

#include <memory> // unique_ptr

#include "MathFunctions.h"
#include "Structure/Database.h"
#include "Structure/SurfaceMesh/MVoronoiCell.h"

// Forward declarations
class Membrane;

/******************************************************************************

The vertex class extends bead,
but it is exclusively used in 2D surface meshwork, and contains
information of its neighbors.

******************************************************************************/

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
        : coord(v), force{}, topoIndex_(topoIndex), mVertex(0)
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

    MVoronoiCell mVertex; // vertex mechanical information

};

#endif
