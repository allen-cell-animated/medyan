#ifndef MEDYAN_Membrane_h
#define MEDYAN_Membrane_h

#include <array>
#include <tuple>
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

    short _memType; // Membrane type

    SubSystem* _subSystem; // SubSystem pointer

    static Database<Membrane*> _membranes; // Collection in SubSystem
    int _Id; // Unique integer id of this membrane

public:
    // Constructors
    // This constructor creates a membrane according to vertex and neighbor data
    Membrane(SubSystem* s, short membraneType,
        vector<tuple<array<double, 3>, vector<size_t>>>& membraneData);

    // This destructor is called when a membrane is to be removed from the system.
    // Removes all vertices, edges and triangles associated with the membrane.
    ~Membrane();

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

    // Print self information
    virtual void printSelf();


};








#endif
