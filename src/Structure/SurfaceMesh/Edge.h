#ifndef MEDYAN_Edge_h
#define MEDYAN_Edge_h

/*
 
An unordered edge contains 2 vertices.
 
*/

#include <array>

#include "common.h"

#include "Trackable.h"
#include "DynamicNeighbor.h"
#include "Component.h"
#include "MathFunctions.h"
#include "Movable.h"

#include "Vertex.h"

// Forward declarations
class Vertex;
class Triangle;
class Compartment;

class Edge:
    public Component,
    public Movable,
    public Trackable {

private:

    size_t _topoIndex; // Index in the meshwork topology.

    static Database<Edge*> _edges; // Collection of edges in SubSystem
    int _id; // Unique integer id of this edge

    void updateCoordinate(); // helper function to update coordiante of this edge

    Compartment* _compartment = nullptr; // The compartment containing this edge

public:
    Edge(Composite *parent, size_t topoIndex);
    ~Edge();

    void setTopoIndex(size_t index) { _topoIndex = index; }

    mathfunc::Vec3 coordinate; // Coordinate of the mid point, updated with updateCoordiante()

    /// Get all instances of this class from the SubSystem
    static const vector<Edge*>& getEdges() {
        return _edges.getElements();
    }
    /// Get Id
    int getId()const { return _id; }

    //@{
    /// Implements Component
    virtual int getType() override { return getParent()->getType(); }
    virtual void printSelf()const override;
    //@}

    //@{
    /// SubSystem management, inherited from Trackable
    virtual void addToSubSystem()override { _edges.addElement(this); }
    virtual void removeFromSubSystem()override { _edges.removeElement(this); }
    //@}

    //@{
    /// Implements Movable
    virtual void updatePosition() override;
    //@}
    Compartment* getCompartment() { return _compartment; }

};

#endif
