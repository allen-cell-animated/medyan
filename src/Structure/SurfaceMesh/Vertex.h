#ifndef MEDYAN_Vertex_h
#define MEDYAN_Vertex_h

#include<vector>

#include "common.h"

#include "Trackable.h"
#include "Movable.h"
#include "Reactable.h"
#include "DynamicNeighbor.h"
#include "Component.h"

#include "Bead.h"

#include "Edge.h"
#include "Triangle.h"
#include "MVoronoiCell.h"

// Forward declarations
class Edge;
class Triangle;

/******************************************************************************

The vertex class extends bead,
but it is exclusively used in 2D surface meshwork, and contains
information of its neighbors.

******************************************************************************/

class Vertex:
    public DynamicNeighbor,
    public Bead {

private:
    unique_ptr<MVoronoiCell> _mVoronoiCell; // pointer to Voronoi cell mechanical information

    // The following vectors on neighbor vertices, triangles, edges must have the same size
    std::vector<Vertex*> _neighborVertices; // tethered neighbors in counter-clockwise direction
    std::vector<Triangle*> _neighborTriangles; // triangles around the vertex. each index i corresponds to
                                       // the triangle (this, _neighborVertices[i], _neighborVertices[i+1])
    std::vector<size_t> _triangleHead; // The index of triangle 0th vertex // TODO: hashmaps
                                       // ...[2] = 1 means the 2nd triangle has vertices at (2, 3, center)
                                       // ...[4] = 2 means the 4th triangle has vertices at (5, center, 4)
    std::vector<Edge*> _neighborEdges; // The tethers, with the same order as neighbor vertices.
    std::vector<size_t> _edgeHead; // The index of edge head vertex
                                   // ...[2] = 0 means the 2nd edge has 0th vertex at center
    
    static Database<Vertex*> _vertices; // Collection of vertices in SubSystem

public:
    ///Main constructor
    Vertex(vector<double> v, Composite* parent, int position);
    
    ///Default constructor
    Vertex(Composite* parent, int position);

    // Get mech Voronoi cell
    MVoronoiCell* getMVoronoiCell() { return _mVoronoiCell.get(); }

    // Get number of tethered neighbors
    size_t getNeighborNum() { return _neighborVertices.size(); }

    // Get tethered neighbors
    std::vector<Vertex*>& getNeighborVertices() { return _neighborVertices; }
    std::vector<Edge*>& getNeighborEdges() { return _neighborEdges; }
    std::vector<Triangle*>& getNeighborTriangles() { return _neighborTriangles; }
    std::vector<size_t>& getEdgeHead() { return _edgeHead; }
    std::vector<size_t>& getTriangleHead() { return _triangleHead; }

    /// Get all instances of this class from the SubSystem
    static const vector<Vertex*>& getVertices() {
        return _vertices.getElements();
    }

};


#endif
