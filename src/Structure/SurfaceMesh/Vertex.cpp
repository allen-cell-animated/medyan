#include "Structure/SurfaceMesh/Vertex.hpp"

#include "Structure/SurfaceMesh/Membrane.hpp" // Membrane::getMesh()

Vertex::Vertex(const Bead::coordinate_type& v, Composite* parent, size_t topoIndex):
    Bead(mathfunc::vec2Vector(v), parent, 0), _topoIndex(topoIndex) {
    
#ifdef MECHANICS
    // eqArea cannot be obtained at this moment
    _mVoronoiCell = std::make_unique<MVoronoiCell>(getType());
#endif

}
