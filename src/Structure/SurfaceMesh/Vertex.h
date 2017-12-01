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
#include "MVoronoiCell.h"
#include "Triangle.h"

/*

The vertex class is just a container that wraps a bead,
but it is exclusively used in 2D surface meshwork, and contains
information of its neighbors.

*/

class Vertex:
    public Component,
    public Trackable,
    public Movable,
    public DynamicNeighbor {

private:
    // Pointers to the bead.
    Bead* _b;

    unique_ptr<MVoronoiCell> _mVoronoiCell; // pointer to Voronoi cell mechanical information

    // The following vectors on neighbor vertices, triangles, edges must have the same size
    std::vector<Vertex*> _tethered; // tethered neighbors in counter-clockwise direction
    std::vector<Triangle*> _triangles; // triangles around the vertex. each index i corresponds to
                                       // the triangle (this, _tethered[i], _tethered[i+1])
    std::vector<Edge*> _edges; // The tethers, with the same order as neighbor vertices.

public:
    Vertex(Composite *parent, Bead *b);

    // Get Beads
    Bead* getBead() { return _b; }
    
    // Get mech Voronoi cell
    MVoronoiCell* getMVoronoiCell() { return _mVoronoiCell.get(); }

    // Get number of tethered neighbors
    size_t getNeighborNum() { return _tethered.size(); }

    // Get tethered neighbor vertices
    std::vector<Vertex*>& getTetheredNeighbors() { return _tethered; }
};


#endif
