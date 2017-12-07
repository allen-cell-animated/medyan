#ifndef MEDYAN_Membrane_h
#define MEDYAN_Membrane_h

#include <vector>

#include "Database.h"
#include "Trackable.h"
#include "Composite.h"

// FORWARD DECLARATIONS
class SubSystem;
class Triangle;
class Edge;
class Vertex;

class Membrane: public Composite, public Trackable {

    friend class Controller;

private:

    vector<Triangle*> _triangleVector; // collection of triangles
    vector<Edge*> _edgeVector; // collection of edges
    vector<Vertex*> _vertexVector; // collection of vertices

    SubSystem* _subSystem; // SubSystem pointer

    static Database<Membrane*> _membranes; // Collection in SubSystem

public:
    // TODO: constructors

    /// Get vector of triangles/edges/vertices that this membrane contains.
    vector<Triangle*>& getTriangleVector() {return _triangleVector;}
    vector<Edge*>& getEdgeVector() { return _edgeVector; }
    vector<Vertex*>& getVertexVector() { return _vertexVector; }

    // SubSystem management, inherited from Trackable
    virtual void addToSubSystem() { _membranes.addElement(this); }
    virtual void removeFromSubSystem() { _membranes.removeElement(this); }
    
    /// Get all instances of this class from the SubSystem
    static const vector<Membrane*>& getMembranes() {
        return _membranes.getElements();
    }
    /// Get the number of membranes in this system
    static int numFilaments() {
        return _membranes.countElements();
    }


};








#endif
