#include "Vertex.h"

Database<Vertex*> Vertex::_vertices;

Vertex::Vertex(vector<double> v, Composite* parent, size_t numNeighbors):
    Bead(v, parent, 0), _id(_vertices.getID()),
    _neighborVertices(numNeighbors, nullptr),
    _neighborTriangles(numNeighbors, nullptr), _triangleHead(numNeighbors),
    _neighborEdges(numNeighbors, nullptr), _edgeHead(numNeighbors) {
    
    _gVoronoiCell = unique_ptr<GVoronoiCell>(new GVoronoiCell(numNeighbors));
    _gVoronoiCell->setVertex(this);
#ifdef MECHANICS
    // eqArea cannot be obtained at this moment
    _mVoronoiCell = unique_ptr<MVoronoiCell>(new MVoronoiCell(getType()));
    _mVoronoiCell->setVertex(this);
#endif

}

void Vertex::updateGeometry(bool calcDerivative, double d) {

    _gVoronoiCell->updateGeometry(calcDerivative, d);

}
