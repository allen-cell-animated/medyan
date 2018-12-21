#include "Vertex.h"

Database<Vertex*> Vertex::_vertices;

Vertex::Vertex(vector<double> v, Composite* parent, size_t topoIndex):
    Bead(v, parent, 0), _id(_vertices.getID()), _topoIndex(topoIndex) {
    
    _gVoronoiCell = unique_ptr<GVoronoiCell>(new GVoronoiCell(getNeighborNum()));
    _gVoronoiCell->setVertex(this);
#ifdef MECHANICS
    // eqArea cannot be obtained at this moment
    _mVoronoiCell = unique_ptr<MVoronoiCell>(new MVoronoiCell(getType()));
    _mVoronoiCell->setVertex(this);
#endif

}
