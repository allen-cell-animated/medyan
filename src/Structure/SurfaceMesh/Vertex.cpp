#include "Structure/SurfaceMesh/Vertex.h"

#include "Structure/SurfaceMesh/Membrane.hpp" // Membrane::getMesh()

Database<Vertex*> Vertex::_vertices;

Vertex::Vertex(vector<double> v, Composite* parent, size_t topoIndex):
    Bead(v, parent, 0), _topoIndex(topoIndex) {
    
#ifdef MECHANICS
    // eqArea cannot be obtained at this moment
    _mVoronoiCell = std::make_unique<MVoronoiCell>(getType());
#endif

}
