#ifndef MEDYAN_Vertex_h
#define MEDYAN_Vertex_h

#include<unordered_map>
#include<vector>

#include "common.h"

#include "Bead.h"

#include "Edge.h"
#include "Triangle.h"
#include "GVoronoiCell.h"
#include "MVoronoiCell.h"

// Forward declarations
class Edge;
class Triangle;
class Membrane;

/******************************************************************************

The vertex class extends bead,
but it is exclusively used in 2D surface meshwork, and contains
information of its neighbors.

******************************************************************************/

class Vertex:
    public Bead // Inherited from bead to receive full features like coordinate and forces.
                // But note that in terms of tracking, when the vertex is added to the system,
                // the base class Bead should also be added to its own collection,
                // i.e. both the bead and the vertex collection should both have the collection.
                // So when initialized, this class
    {
    
    friend class Membrane; // Membrane class can manage id of this vertex

private:
    size_t _topoIndex; // Index in the meshwork topology.

    unique_ptr<GVoronoiCell> _gVoronoiCell; // pointer to Voronoi cell geometric information
    unique_ptr<MVoronoiCell> _mVoronoiCell; // pointer to Voronoi cell mechanical information

    // The following vectors on neighbor vertices, triangles, edges must have the same size
    std::vector<Vertex*> _neighborVertices; // tethered neighbors in counter-clockwise direction
    std::unordered_map<Vertex*, size_t> _neighborVertexIndices; // maps vertex pointers back to indices
    std::vector<Triangle*> _neighborTriangles; // triangles around the vertex. each index i corresponds to
                                       // the triangle (this, _neighborVertices[i], _neighborVertices[i+1])
    std::vector<size_t> _triangleHead; // The index of triangle 0th vertex
                                       // ...[2] = 1 means the 2nd triangle has vertices at (2, 3, center)
                                       // ...[4] = 2 means the 4th triangle has vertices at (5, center, 4)
    std::vector<Edge*> _neighborEdges; // The tethers, with the same order as neighbor vertices.
    std::vector<size_t> _edgeHead; // The index of edge head vertex
                                   // ...[2] = 0 means the 2nd edge has 0th vertex at center
    
    static Database<Vertex*> _vertices; // Collection of vertices in SubSystem
    int _id; // Unique integer id of this vertex
             // Warning: By 2017/12/18 the Bead::_ID attribute is never initialized or used,
             // but if the beads need to be assigned an ID at construction, one should
             // know that the vertex id either is DIFFERENT from or HIDES the bead id,
             // so that the bead id is either NOT USED or HIDDEN by the id in vertex structure.
    
    size_t _membraneVertexIdx; // The index of this in the _vertexVector of the parent membrane

public:
    ///Main constructor
    Vertex(vector<double> v, Composite* parent, size_t topoIndex);

    void setTopoIndex(size_t index) { _topoIndex = index; }
    
    // Get geo Voronoi cell
    GVoronoiCell* getGVoronoiCell() { return _gVoronoiCell.get(); }
    // Get mech Voronoi cell
    MVoronoiCell* getMVoronoiCell() { return _mVoronoiCell.get(); }

    // Get number of tethered neighbors
    size_t getNeighborNum()const { return _neighborVertices.size(); }

    // Get tethered neighbors
    std::vector<Vertex*>& getNeighborVertices() { return _neighborVertices; }
    std::unordered_map<Vertex*, size_t>& getNeighborVertexIndices() { return _neighborVertexIndices; }
    std::vector<Edge*>& getNeighborEdges() { return _neighborEdges; }
    std::vector<Triangle*>& getNeighborTriangles() { return _neighborTriangles; }
    std::vector<size_t>& getEdgeHead() { return _edgeHead; }
    std::vector<size_t>& getTriangleHead() { return _triangleHead; }

    /// Get all instances of this class from the SubSystem
    static const vector<Vertex*>& getVertices() {
        return _vertices.getElements();
    }
    /// Get ID
    int getId()const { return _id; }
    size_t getMembraneVertexIdx()const { return _membraneVertexIdx; }

    //@{
    /// SubSystem management, inherited from Bead: Trackable
    // Incremental to the functions in Bead class.
    virtual void addToSubSystem()override {
        _vertices.addElement(this);
        Bead::addToSubSystem();
    }
    virtual void removeFromSubSystem()override {
        _vertices.removeElement(this);
        Bead::removeFromSubSystem();
    }
    //@}


};


#endif
