#include "Membrane.h"

#include "common.h"

#include "SubSystem.h"

#include "Triangle.h"
#include "Edge.h"
#include "Vertex.h"

#include "MathFunctions.h"
using namespace mathfunc;

Database<Membrane*> Membrane::_membranes;

Membrane::Membrane(SubSystem* s, short membraneType,
    vector<tuple<array<double, 3>, vector<size_t>>>& membraneData):
    Trackable(), _subSystem(s), _memType(membraneType), _Id(_membranes.getID()) {
    
    // Build the meshwork using vertex and neighbor information

    size_t numVertices = membraneData.size();
    if(numVertices == 0) return;

    /**************************************************************************
        Setting up vertices and neighbors
    **************************************************************************/
    // Add the vertices
    for(auto& vertexData : membraneData) {
        _subSystem->addTrackable<Vertex>(
            array2Vector<double, 3>(get<0>(vertexData)),
            this,
            get<1>(vertexData).size()
        );
    }
    int futureIdx = Vertex::_vertices.getElements().size();
    int firstIdx = futureIdx - numVertices;

    // Register the neighbors
    for(int idx = firstIdx; idx < futureIdx; ++idx) {
        auto& neighborData = get<1>(membraneData[idx - firstIdx]);
        Vertex* centerVertex = Vertex::_vertices.getElements()[idx];
        size_t numNeighbors = centerVertex->getNeighborNum();
        for(size_t nIdx = 0; nIdx < numNeighbors; ++nIdx) {
            centerVertex->getNeighborVertices()[nIdx] = Vertex::_vertices.getElements()[neighborData[nIdx] + firstIdx];
        }
    }

    /**************************************************************************
        Setting up edges and triangles
    **************************************************************************/
    for(int idx = firstIdx; idx < futureIdx; ++idx) {
        Vertex* centerVertex = Vertex::_vertices.getElements()[idx];
        size_t numNeighbors = centerVertex->getNeighborNum();
        for(int nIdx = 0; nIdx < numNeighbors; ++nIdx) {

            Vertex* nVertex = centerVertex->getNeighborVertices()[nIdx];
            Vertex* nnVertex = centerVertex->getNeighborVertices()[(nIdx + 1) % numNeighbors];

            // Edge registration
            if(centerVertex->getNeighborEdges()[nIdx] == nullptr) { // Edge not registered
                _subSystem->addTrackable<Edge>(
                    this,
                    centerVertex,
                    centerVertex->getNeighborVertices()[nIdx]
                );
                Edge* lastAddedEdge = Edge::_edges.getElements().back();
                centerVertex->getNeighborEdges()[nIdx] = lastAddedEdge;
                nVertex->getNeighborEdges()[nIdx]; // TODO: Neighbor structure
            }
        }
    }
    // TODO: Implement this

}

Membrane::~Membrane() {
    for(auto& v: _vertexVector) _subSystem->removeTrackable<Vertex>(v);
    for(auto& e: _edgeVector) _subSystem->removeTrackable<Edge>(e);
    for(auto& t: _triangleVector) _subSystem->removeTrackable<Triangle>(t);
}

void Membrane::printSelf() {
    
    cout << endl;
    
    cout << "Membrane: ptr = " << this << endl;
    cout << "Membrane Id = " << _Id << endl;
    cout << "Membrane type = " << _memType << endl;
    
    cout << endl;
    cout << "Triangle information..." << endl;
    
    for(auto t : _triangleVector)
        t->printSelf();
    
    cout << endl;
    
}
