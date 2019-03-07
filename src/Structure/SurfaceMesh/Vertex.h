#ifndef MEDYAN_Vertex_h
#define MEDYAN_Vertex_h

#include<vector>

#include "common.h"

#include "Bead.h"

#include "Edge.h"
#include "Triangle.h"
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

    unique_ptr<MVoronoiCell> _mVoronoiCell; // pointer to Voronoi cell mechanical information
    
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
    
    // Get mech Voronoi cell
    MVoronoiCell* getMVoronoiCell() { return _mVoronoiCell.get(); }

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
