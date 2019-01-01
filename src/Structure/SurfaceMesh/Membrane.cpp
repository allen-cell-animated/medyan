#include "Membrane.hpp"

#include <limits>
#include <stdexcept>
#include <unordered_set>

#include "common.h"
#include "MathFunctions.h"

#include "Compartment.h"
#include "core/controller/GController.h"
#include "SubSystem.h"
#include "SysParams.h"

#include "Triangle.h"
#include "Edge.h"
#include "Vertex.h"
#include "GTriangle.h"
#include "MTriangle.h"
#include "GEdge.h"
#include "GVoronoiCell.h"
#include "MVoronoiCell.h"
#include "GMembrane.h"
#include "MMembrane.h"

using namespace mathfunc;

Database<Membrane*> Membrane::_membranes;

Membrane::Membrane(
    SubSystem* s,
    short membraneType,
    const std::vector< MembraneMeshAttribute::coordinate_type >& vertexCoordinateList,
    const std::vector< std::array< size_t, 3 > >& triangleVertexIndexList
) : Trackable(false, false, false, false, true), Geometric(),
    _mesh(MembraneMeshAttribute::MetaAttribute{s, this}),
    _subSystem(s), _memType(membraneType), _id(_membranes.getID()) {
    
    // Build the meshwork using vertex and triangle information
    _mesh.init(vertexCoordinateList, triangleVertexIndexList);

    size_t numVertices = membraneData.size();
    if(numVertices == 0) return;

    /**************************************************************************
        Setting up vertices and neighbors
    **************************************************************************/
    // Add the vertices
    size_t vertexIndex = 0;
    _vertexVector.reserve(numVertices);
    for(auto& vertexData : membraneData) {
        Vertex* lastAddedVertex = _subSystem->addTrackable<Vertex>(
            array2Vector<double, 3>(get<0>(vertexData)),
            this,
            get<1>(vertexData).size()
        );
        _vertexVector.push_back(lastAddedVertex); // Add to its own storage
        lastAddedVertex->_membraneVertexIdx = vertexIndex++;
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
                Edge* lastAddedEdge = _subSystem->addTrackable<Edge>(this, centerVertex, nVertex);
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

                // Calculate the length of the edge
                lastAddedEdge->getGEdge()->calcLength();
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
                Triangle* lastAddedTriangle = _subSystem->addTrackable<Triangle>(this, centerVertex, nVertex, nnVertex);
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
                lastAddedTriangle->getEdges() = {{
                    centerVertex->getNeighborEdges()[nIdx],
                    nVertex->getNeighborEdges()[idx12],
                    nnVertex->getNeighborEdges()[idx20]
                }};
                lastAddedTriangle->getEdgeHead() = {{
                    centerVertex->getEdgeHead()[nIdx],
                    nVertex->getEdgeHead()[idx12],
                    nnVertex->getEdgeHead()[idx20]
                }};
                
                // Bind the triangle to edges
                for(size_t eIdx = 0; eIdx < 3; ++eIdx)
                    lastAddedTriangle->getEdges()[eIdx]->getTriangles()[lastAddedTriangle->getEdgeHead()[eIdx]] = lastAddedTriangle;

                // Calculate the area of the triangle and set it as eqArea
                lastAddedTriangle->getGTriangle()->calcArea();
#ifdef MECHANICS
                lastAddedTriangle->getMTriangle()->setEqArea(
                    lastAddedTriangle->getGTriangle()->getArea() *
                    SysParams::Mechanics().MemEqAreaFactor[membraneType]
                );
#endif
				// Calculate angles for the use of Voronoi cells
				lastAddedTriangle->getGTriangle()->calcTheta();
            }
        }
    }

    /**************************************************************************
        Setting up Voronoi cells
    **************************************************************************/
    for(size_t idx = 0; idx < numVertices; ++idx) {
        GVoronoiCell* gvc = _vertexVector[idx]->getGVoronoiCell();
        gvc->calcArea();
#ifdef MECHANICS
        MVoronoiCell* mvc = _vertexVector[idx]->getMVoronoiCell();
        // Set the current area as eqArea
        mvc->setEqArea(gvc->getArea() * SysParams::Mechanics().MemEqAreaFactor[membraneType]);
#endif
    }

    /**************************************************************************
        Setting up MMembrane object and find volume
    **************************************************************************/
    _gMembrane = unique_ptr<GMembrane>(new GMembrane);
    _gMembrane->setMembrane(this);
    _gMembrane->calcVolume();
#ifdef MECHANICS
    _mMembrane = unique_ptr<MMembrane>(new MMembrane);
    _mMembrane->setMembrane(this);
    _mMembrane->setEqVolume(_gMembrane->getVolume());
#endif

}

void Membrane::printSelf()const {
    
    cout << endl;
    
    cout << "Membrane: ptr = " << this << endl;
    cout << "Membrane Id = " << _id << endl;
    cout << "Membrane type = " << _memType << endl;
    
    cout << endl;
    
    cout << endl;
    
}

double Membrane::meshworkQuality()const {
    /*
    This function calculates the quality of the meshwork of this membrane, and
    the result is represented as a value between 0 and 1, 0 being worst and 1
    being the best case (all equilateral).

    The criterion used is the minimum angle in each triangle, parametrized and
    averaged over all triangles. And the calculation requires the result of
        - The angle calculation of all triangles
    */

    double res = 0;

    // TODO implementation

    return res;
}
