#include "Structure/SurfaceMesh/Vertex.h"

#include "Structure/SurfaceMesh/Membrane.hpp" // Membrane::getMesh()

Database<Vertex*> Vertex::_vertices;

Vertex::Vertex(vector<double> v, Composite* parent, size_t topoIndex):
    Bead(v, parent, 0), _id(_vertices.getID()), _topoIndex(topoIndex) {
    
#ifdef MECHANICS
    // eqArea cannot be obtained at this moment
    _mVoronoiCell = std::make_unique<MVoronoiCell>(getType());
#endif

}

size_t Vertex::getNeighborNum()const {
    const auto& mesh = static_cast<const Membrane*>(getParent())->getMesh();
    size_t num = 0;
    mesh.forEachHalfEdgeTargetingVertex(_topoIndex, [&num](size_t hei) { ++num; });
    return num;
}
