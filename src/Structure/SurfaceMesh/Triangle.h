#ifndef MEDYAN_Triangle_h
#define MEDYAN_Triangle_h

#include <array>

#include "common.h"

#include "Database.h"
#include "Trackable.h"
#include "Movable.h"
//#include "Reactable.h"
#include "DynamicNeighbor.h"
#include "Component.h"

#include "Vertex.h"
#include "Edge.h"
#include "MTriangle.h"

// Forward declarations
class Edge;
class Vertex;
class Compartment;

/******************************************************************************
Triangles are the only element in the meshwork that has area and act as patches
that construct the surface.

The triangle patches have geometric and mechanical properties.

The Triangle class has pointers to the vertices and edges.
******************************************************************************/
class Triangle:
    public Component,
    public Trackable,
    public Movable,
    // public Reactable,
    public DynamicNeighbor {

private:
    size_t _topoIndex; // Index in the meshwork topology.

    unique_ptr<MTriangle> _mTriangle; // pointer to mech triangle

    static Database<Triangle*> _triangles; // Collection of triangles in SubSystem
    int _id; // Unique integer id of this triangle

    void updateCoordinate(); // helper function to update coordiante of this triangle

    Compartment* _compartment = nullptr; // The compartment containing this triangle

public:
    std::array<double, 3> coordinate; // Coordinate of the center point, updated with updateCoordiante()

    Triangle(Composite *parent, size_t topoIndex);

    void setTopoIndex(size_t index) { _topoIndex = index; }
    size_t getTopoIndex() const { return _topoIndex; }

    // Get mech triangle
    MTriangle* getMTriangle() { return _mTriangle.get(); }

    /// Get all instances of this class from the SubSystem
    static const vector<Triangle*>& getTriangles() {
        return _triangles.getElements();
    }
    /// Get ID
    int getId()const { return _id; }

    //@{
    /// SubSystem management, inherited from Trackable
    virtual void addToSubSystem()override { _triangles.addElement(this); }
    virtual void removeFromSubSystem()override { _triangles.removeElement(this); }
    //@}

    //@{
    /// Implements Component
    virtual int getType() override { return getParent()->getType(); }
    virtual void printSelf()const override;
    //@}

    //@{
    /// Implements Movable
    virtual void updatePosition() override;
    //@}
    Compartment* getCompartment() { return _compartment; }


};

#endif
