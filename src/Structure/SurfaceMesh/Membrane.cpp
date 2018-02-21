#include "Membrane.h"

#include <limits>
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
#include "GMembrane.h"
#include "MMembrane.h"

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
    if(calcDerivative) _gMembrane->calcVolume(); else _gMembrane->calcStretchedVolume(d);
}

double Membrane::signedDistance(const std::array<double, 3>& p, bool safe)const {
    if(!_isClosed) throw std::logic_error("Membrane is not closed while trying to find signed distance field.");

    /**************************************************************************
    The function works in the following procedure:

    - Iterate through a certain amount of triangles, and for each triangle
        - Find the projection of the point on the triangle plane, and determine
          which element is responsible for being the closest to the point
          (which vertex/ which edge/ this triangle).
        - Find the unsigned distance with the closest element, and then find
          the signed distance using the normal or pseudo normal. Record the
          value with the smallest unsigned distance.
    
    Before this function is used, the following must be calculated:
        - The positions of all the elements are updated
        - The normal and pseudo normal at the triangles, edges and vertices
        - The length of edges
    
    In fact, the signed distance field serves as a good candidate for membrane
    boundary potential. However, this field is only C0-continuous, so some
    configurations might result in non-stop situation in CG method.
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
        vector<int> indices = GController::getCompartmentIndices(mathfunc::array2Vector(p));

        if(indices.size() != 3)
            throw std::logic_error("Only 3D compartments are allowed in the membrane signed distance calculation.");

		unordered_set<Triangle*> triSet; // Set of triangles to be searched

        // Find all the neighbor compartments and find triangles to be added to search set.
        // The for loops are written like this for aesthetic reasons. Be very careful with the scope below.
        for(int v0 = -1; v0 <= 1; ++v0) for(int v1 = -1; v1 <= 1; ++v1) for(int v2 = -1; v2 <= 1; ++v2) {

            vector<int> newIndices = {indices[0] + v0, indices[1] + v1, indices[2] + v2};
            if(GController::indicesOutOfBound(indices)) {
                // Compartment not found. Simply ignore it.
                continue;
            }
            
            Compartment* c = GController::getCompartment(newIndices);

            // Compartment exists
            const unordered_set<Triangle*>& triSetEach = c->getTriangles();
            triSet.insert(triSetEach.begin(), triSetEach.end());
        }

        if(triSet.empty()) {
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

    double minAbsDistance = numeric_limits<double>::infinity();
    for(Triangle* t: *loopingTriangles) {
        /**********************************************************************
        Calculate the barycentric coordinate of the projection point p'

        See Heidrich 2005, Computing the Barycentric Coordinates of a Projected
        Point.
        **********************************************************************/
        const array<double, 3> v0 = vector2Array<double, 3>(t->getVertices()[0]->coordinate);
        const array<double, 3> v1 = vector2Array<double, 3>(t->getVertices()[1]->coordinate);
        const array<double, 3> v2 = vector2Array<double, 3>(t->getVertices()[2]->coordinate);

        array<double, 3> n = vectorProduct(v0, v1, v0, v2);
        double oneOver4AreaSquared = 1.0 / dotProduct(n, n);

        double b1 = dotProduct(vectorProduct(v0, p, v0, v2), n) * oneOver4AreaSquared;
        double b2 = dotProduct(vectorProduct(v0, v1, v0, p), n) * oneOver4AreaSquared;
        double b0 = 1.0 - b1 - b2;

        // Now p' = b0*v0 + b1*v1 + b2*v2

        double d = numeric_limits<double>::infinity();
        // b0, b1 and b2 cannot all be negative at the same time, leaving 7 possible combinations
        if(b0 >= 0) {
            if(b1 >= 0) {
                if(b2 >= 0) {
                    // Inside triangle
                    d = dotProduct(t->getGTriangle()->getUnitNormal(), vectorDifference(p, v0));
                }
                else {
                    // On edge 01
                    auto ge = t->getEdges()[0]->getGEdge();
                    d = magnitude(vectorProduct(v0, p, v0, v1)) / ge->getLength();
                    if(dotProduct(ge->getPseudoUnitNormal(), vectorDifference(p, v0)) < 0) d = -d;
                }
            }
            else {
                if(b2 >= 0) {
                    // On edge 02
                    auto ge = t->getEdges()[2]->getGEdge();
                    d = magnitude(vectorProduct(v0, p, v0, v2)) / ge->getLength();
                    if(dotProduct(ge->getPseudoUnitNormal(), vectorDifference(p, v0)) < 0) d = -d;
                }
                else {
                    // On vertex 0
                    d = twoPointDistance(v0, p);
                    if(dotProduct(t->getVertices()[0]->getGVoronoiCell()->getPseudoUnitNormal(), vectorDifference(p, v0)) < 0)
                        d = -d;
                }
            }
        }
        else {
            if(b1 >= 0) {
                if(b2 >= 0) {
                    // On edge 12
                    auto ge = t->getEdges()[1]->getGEdge();
                    d = magnitude(vectorProduct(v1, p, v1, v2)) / ge->getLength();
                    if(dotProduct(ge->getPseudoUnitNormal(), vectorDifference(p, v1)) < 0) d = -d;
                }
                else {
                    // On vertex 1
                    d = twoPointDistance(v1, p);
                    if(dotProduct(t->getVertices()[1]->getGVoronoiCell()->getPseudoUnitNormal(), vectorDifference(p, v1)) < 0)
                        d = -d;
                }
            }
            else {
                if(b2 >= 0) {
                    // On vertex 2
                    d = twoPointDistance(v2, p);
                    if(dotProduct(t->getVertices()[2]->getGVoronoiCell()->getPseudoUnitNormal(), vectorDifference(p, v2)) < 0)
                        d = -d;
                }
                else {
                    // The program should never go here
                    throw logic_error("Hey, you know the barycentric coordinates cannot be all negative. Something has gone terribly wrong.");
                }
            }
        }

        // Update with distance with less absolute value
        if(abs(d) < abs(minAbsDistance)) minAbsDistance = d;
    }
    
    return minAbsDistance;
}

bool Membrane::contains(const std::array<double, 3>& point)const {
    return signedDistance(point, true) < 0.0;
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
