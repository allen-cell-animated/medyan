#include "Membrane.h"

#include <stdexcept>
#include <unordered_set>

#include "common.h"
#include "MathFunctions.h"
using namespace mathfunc;

#include "SubSystem.h"
#include "Compartment.h"
#include "GController.h"

#include "Triangle.h"
#include "Edge.h"
#include "Vertex.h"
#include "GTriangle.h"
#include "MTriangle.h"
#include "GEdge.h"
#include "GVoronoiCell.h"
#include "MVoronoiCell.h"

Database<Membrane*> Membrane::_membranes;

Membrane::Membrane(SubSystem* s, short membraneType,
    vector<tuple<array<double, 3>, vector<size_t>>>& membraneData):
    Trackable(false, false, false, false, true), Geometric(),
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
                
                // Bind the triangle to edges
                for(size_t eIdx = 0; eIdx < 3; ++eIdx)
                    lastAddedTriangle->getEdges()[eIdx]->getTriangles()[lastAddedTriangle->getEdgeHead()[eIdx]] = lastAddedTriangle;

                // Calculate the area of the triangle and set it as eqArea
                lastAddedTriangle->getGTriangle()->calcArea();
#ifdef MECHANICS
                lastAddedTriangle->getMTriangle()->setEqArea(lastAddedTriangle->getGTriangle()->getArea());
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
        mvc->setEqArea(gvc->getArea());
#endif
    }

}

Membrane::~Membrane() {
    for(auto& v: _vertexVector) _subSystem->removeTrackable<Vertex>(v);
    for(auto& e: _edgeVector) _subSystem->removeTrackable<Edge>(e);
    for(auto& t: _triangleVector) _subSystem->removeTrackable<Triangle>(t);
}

void Membrane::printSelf()const {
    
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
    /**************************************************************************
    Updates the geometric properties of all the elements of the membrane. This
    MUST be called before any energy calculation of the membrane.

    Due to the geometry dependencies, the operation order below is important.

    Some of the geometric properties are not needed in energy computation, so
    for efficiency, they are not updated in this function. They must be updated
    manually before use. They are:
        - <None>
    **************************************************************************/

    for(auto& e: _edgeVector) {
        auto ge = e->getGEdge();
        if(calcDerivative) ge->calcLength(); else ge->calcStretchedLength(d);
    }
    for(auto& t: _triangleVector) {
        auto gt = t->getGTriangle();
        if(calcDerivative) gt->calcTheta(); else gt->calcStretchedTheta(d);
        if(calcDerivative) gt->calcArea(); else gt->calcStretchedArea(d);
        if(calcDerivative) gt->calcUnitNormal(); else gt->calcStretchedUnitNormal(d);
    }
    for(auto& e: _edgeVector) {
        auto ge = e->getGEdge();
        if(calcDerivative) ge->calcPseudoUnitNormal(); else ge->calcStretchedPseudoUnitNormal(d);
    }
    for(auto& v: _vertexVector) {
        auto gv = v->getGVoronoiCell();
        if(calcDerivative) gv->calcArea(); else gv->calcStretchedArea(d);
        if(calcDerivative) gv->calcCurv(); else gv->calcStretchedCurv(d);
        if(calcDerivative) gv->calcPseudoUnitNormal(); else gv->calcStretchedPseudoUnitNormal(d);
    }
}

double Membrane::signedDistance(const std::array<double, 3>& p, bool safe)const {
    if(!_isClosed) throw std::logic_error("Membrane is not closed while trying to find signed distance field.");

    /**************************************************************************
    The functions works in the following procedure.

    - Iterate through a certain amount of triangles, and for each triangle
        - Find the projection of the point on the triangle plane, and determine
          which element is responsible for being the closest to the point
          (which vertex/ which edge/ this triangle).
        - Find the unsigned distance with the closest element. Record the
          the value and the element with the smallest such distance found.
    - With the element with the smallest unsigned distance, use the normal or
      pesudonormal to find the signed distance.
    
    Before this function is used, it is required that:
        - The positions of all the elements are updated
    **************************************************************************/

    /**************************************************************************
    Determine the triangles needed to be looped through

    Safe        : Use all triangles
    Not safe    : Search neighboring 27 compartments for triangles
    **************************************************************************/
    const vector<Triangle*>* loopingTriangles = nullptr;

    unique_ptr<vector<Triangle*>> limitedLoopingTriangles;
    if(safe) {
        // Loop through all the triangles
        loopingTriangles = &_triangleVector;
    }
    else {
        // Find the current compartment indices containing point p
        vector<size_t> indices;
        try { indices = GController::getCompartmentIndices(mathfunc::array2Vector(p)); }
        catch (exception& e) {
            cout << e.what() << endl;
            printSelf();
            exit(EXIT_FAILURE);
        }

        if(indices.size() != 3)
            throw std::logic_error("Only 3D compartments are allowed in the membrane signed distance calculation.");

		unordered_set<Triangle*> triSet; // Set of triangles to be searched

        // Find all the neighbor compartments and find triangles to be added to search set.
        // The for loops are written like this for aesthetic reasons. Be very careful with the scope below.
        for(int v0: {-1, 0, 1}) for(int v1: {-1, 0, 1}) for(int v2: {-1, 0, 1}) {

            vector<size_t> newIndices = {indices[0] + v0, indices[1] + v1, indices[2] + v2};
            Compartment* c = nullptr;
            try { c = GController::getCompartment(newIndices); }
            catch(const OutOfBoundsException& e) {
                // Compartment not found. Simply ignore it.
                continue;
            }

            // Compartment exists
            const unordered_set<Triangle*>& triSetEach = c->getTriangles();
            triSet.insert(triSetEach.begin(), triSetEach.end());
        }

        if(triSet.empty) {
            cout << "Warning: triangles not found in neighboring compartments. "
                 << "Getting the result from safe mode."
                 << endl;
            return signedDistance(p, true);
        }

        // Copy the triangles in the unordered set into a vector
        limitedLoopingTriangles = unique_ptr<vector<Triangle*>>(new vector<Triangle*>(triSet.begin(), triSet.end()));

        // Make the new vector the vector to be looped
        loopingTriangles = limitedLoopingTriangles.get();
    }

    /**************************************************************************
     * TODO
    **************************************************************************/
    for(Triangle* t: *loopingTriangles) {

    }
    
    // Helper matrix for finding the projection point P' on the triangle plane with a given point P
    // First find the projection of P on two edges 0--1 and 0--2 as a vector in 2D
    //   b = (dot(r0P, r01), dot(r0P, r02))' = (r01 r02)' * r0P
    //                                         ^~~~~~~~~~
    //                                         2x3 matrix
    // Then apply the _projMat (A) to find vector c: c = A * b
    // Then r0P' = c[0] * r01 + c[1] * r02 = (r01 r02) * c
    //                                       ^~~~~~~~~
    //                                       3x2 matrix
    std::array<std::array<double, 2>, 2> _projMat; // The helper matrix for finding the

    return false; // TODO: Implement it
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

    for(Triangle* t: _triangleVector) {
        auto gt = t->getGTriangle();

        // Find minimum angle
        double minAngle = M_PI / 3; // The largest possible minimum angle in a triangle
        for(double eachTheta: gt->getTheta())
            if(eachTheta < minAngle) minAngle = eachTheta;
        
        // Parametrize the minimum angle and add to result
        double q = minAngle * 3 / M_PI; // Normalize, so that q should be within 0 and 1
        q *= q; // Small angles get more "emphasized"
        res += q;
    }

    res /= _triangleVector.size(); // Average instead of sum

    return res;
}
