#include "Membrane.h"

#include <stdexcept>

#include "common.h"

#include "SubSystem.h"

#include "Triangle.h"
#include "Edge.h"
#include "Vertex.h"
#include "MTriangle.h"
#include "MEdge.h"
#include "MVoronoiCell.h"

#include "MathFunctions.h"
using namespace mathfunc;

Database<Membrane*> Membrane::_membranes;

Membrane::Membrane(SubSystem* s, short membraneType,
    vector<tuple<array<double, 3>, vector<size_t>>>& membraneData):
    Trackable(false, false, false, false, true), Geometric(), // Self management of geometric behavior
    _subSystem(s), _memType(membraneType), _id(_membranes.getID()) {
    
    // Build the meshwork using vertex and neighbor information

    size_t numVertices = membraneData.size();
    if(numVertices == 0) return;

    /**************************************************************************
        Setting up vertices and neighbors
    **************************************************************************/
    // Add the vertices
    _vertexVector.reserve(numVertices);
    for(auto& vertexData : membraneData) {
        _subSystem->addTrackable<Vertex>(
            array2Vector<double, 3>(get<0>(vertexData)),
            this,
            get<1>(vertexData).size()
        );
        _vertexVector.push_back(Vertex::getVertices().back()); // Add to its own storage
    }

    // Register the neighbors
    for(int idx = 0; idx < numVertices; ++idx) {
        auto& neighborData = get<1>(membraneData[idx]);
        Vertex* centerVertex = _vertexVector[idx];
        size_t numNeighbors = centerVertex->getNeighborNum();
        for(size_t nIdx = 0; nIdx < numNeighbors; ++nIdx) {
            Vertex* nVertex = _vertexVector[neighborData[nIdx]];
            centerVertex->getNeighborVertices()[nIdx] = nVertex;
            centerVertex->getNeighborVertexIndices()[nVertex] = nIdx;
        }
    }

    /**************************************************************************
        Setting up edges
    **************************************************************************/
    for(int idx = 0; idx < numVertices; ++idx) {
        Vertex* centerVertex = _vertexVector[idx];
        size_t numNeighbors = centerVertex->getNeighborNum();
        for(int nIdx = 0; nIdx < numNeighbors; ++nIdx) {

            Vertex* nVertex = centerVertex->getNeighborVertices()[nIdx];

            // Edge registration
            if(centerVertex->getNeighborEdges()[nIdx] == nullptr) { // Edge not registered
                _subSystem->addTrackable<Edge>(this, centerVertex, nVertex);
                Edge* lastAddedEdge = Edge::getEdges().back();
                _edgeVector.push_back(lastAddedEdge); // Add to its own storage

                // Bind the edge to vertices, and check whether neighbor exists
                size_t backToCenterIdx = 0;
                try {
                    backToCenterIdx = nVertex->getNeighborVertexIndices().at(centerVertex);
                }
                catch(const std::out_of_range& oor) {
                    cout << "An error occured when trying to add edges of the meshwork. "
                         << "Neighbors must pair with each other. Exiting."
                         << endl;
                    exit(EXIT_FAILURE);
                }

                centerVertex->getNeighborEdges()[nIdx] = lastAddedEdge;
                nVertex->getNeighborEdges()[backToCenterIdx] = lastAddedEdge;
                
                centerVertex->getEdgeHead()[nIdx] = 0;
                nVertex->getEdgeHead()[backToCenterIdx] = 1;

#ifdef MECHANICS
                // Calculate the length of the edge
                lastAddedEdge->getMEdge()->calcLength();
#endif
            }
        }
    }

    /**************************************************************************
        Setting up triangles
    **************************************************************************/
    for(int idx = 0; idx < numVertices; ++idx) {
        Vertex* centerVertex = _vertexVector[idx];
        size_t numNeighbors = centerVertex->getNeighborNum();
        for(int nIdx = 0; nIdx < numNeighbors; ++nIdx) {

            Vertex* nVertex = centerVertex->getNeighborVertices()[nIdx];
            Vertex* nnVertex = centerVertex->getNeighborVertices()[(nIdx + 1) % numNeighbors];

            // Triangle registration
            if(centerVertex->getNeighborTriangles()[nIdx] == nullptr) { // Triangle not registered
                _subSystem->addTrackable<Triangle>(this, centerVertex, nVertex, nnVertex);
                Triangle* lastAddedTriangle = Triangle::getTriangles().back();
                _triangleVector.push_back(lastAddedTriangle); // Add to its own storage

                // Bind the triangle to vertices, and check whether neighbor exists
                size_t idx12 = 0, idx20 = 0;
                try {
                    idx12 = nVertex->getNeighborVertexIndices().at(nnVertex);
                    idx20 = nnVertex->getNeighborVertexIndices().at(centerVertex);
                }
                catch(const std::out_of_range& oor) {
                    cout << "An error occured when trying to add triangles of the meshwork. "
                         << "Vertices should be able to form triangles. Exiting."
                         << endl;
                    exit(EXIT_FAILURE);
                }

                centerVertex->getNeighborTriangles()[nIdx] = lastAddedTriangle;
                nVertex->getNeighborTriangles()[idx12] = lastAddedTriangle;
                nnVertex->getNeighborTriangles()[idx20] = lastAddedTriangle;
                
                centerVertex->getTriangleHead()[nIdx] = 0;
                nVertex->getTriangleHead()[idx12] = 2;
                nnVertex->getTriangleHead()[idx20] = 1;

                // Bind edges to the triangle
                lastAddedTriangle->getEdges() = {
                    centerVertex->getNeighborEdges()[nIdx],
                    nVertex->getNeighborEdges()[idx12],
                    nnVertex->getNeighborEdges()[idx20]
                };
                lastAddedTriangle->getEdgeHead() = {
                    centerVertex->getEdgeHead()[nIdx],
                    nVertex->getEdgeHead()[idx12],
                    nnVertex->getEdgeHead()[idx20]
                };

#ifdef MECHANICS
                // Calculate the area of the triangle and set it as eqArea
                lastAddedTriangle->getMTriangle()->calcArea();
                lastAddedTriangle->getMTriangle()->setEqArea(lastAddedTriangle->getMTriangle()->getArea());
				// Calculate angles for the use of Voronoi cells
				lastAddedTriangle->getMTriangle()->calcTheta();
#endif
            }
        }
    }

#ifdef MECHANICS
    /**************************************************************************
        Setting up Voronoi cells
    **************************************************************************/
    for(int idx = 0; idx < numVertices; ++idx) {
        MVoronoiCell* mvc = _vertexVector[idx]->getMVoronoiCell();
        // Calculate area and set it as eqArea
        mvc->calcArea();
        mvc->setEqArea(mvc->getArea());
    }
#endif

}

Membrane::~Membrane() {
    for(auto& v: _vertexVector) _subSystem->removeTrackable<Vertex>(v);
    for(auto& e: _edgeVector) _subSystem->removeTrackable<Edge>(e);
    for(auto& t: _triangleVector) _subSystem->removeTrackable<Triangle>(t);
}

void Membrane::printSelf() {
    
    cout << endl;
    
    cout << "Membrane: ptr = " << this << endl;
    cout << "Membrane Id = " << _id << endl;
    cout << "Membrane type = " << _memType << endl;
    
    cout << endl;
    cout << "Triangle information..." << endl;
    
    for(auto t : _triangleVector)
        t->printSelf();
    
    cout << endl;
    
}

void Membrane::updateGeometry(bool calcDerivative, double d) {
    for(auto& e: _edgeVector) e->updateGeometry(calcDerivative, d);
    for(auto& t: _triangleVector) t->updateGeometry(calcDerivative, d);
    for(auto& v: _vertexVector) v->updateGeometry(calcDerivative, d);
}
