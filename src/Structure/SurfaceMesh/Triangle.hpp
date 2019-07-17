#ifndef MEDYAN_Structure_SurfaceMesh_Triangle_Hpp
#define MEDYAN_Structure_SurfaceMesh_Triangle_Hpp

#include "common.h"
#include "Database.h"
#include "Trackable.h"
#include "Movable.h"
//#include "Reactable.h"
#include "DynamicNeighbor.h"
#include "MathFunctions.h"
#include "Structure/SurfaceMesh/MTriangle.hpp"

// Forward declarations
class Compartment;
class Membrane;

/******************************************************************************
Triangles are the only element in the meshwork that has area and act as patches
that construct the surface.

The triangle patches have geometric and mechanical properties.

The Triangle class has pointers to the vertices and edges.
******************************************************************************/
class Triangle:
    public Trackable,
    public Movable,
    // public Reactable,
    public DynamicNeighbor,
    public Database< Triangle, false > {

private:
    Membrane* _parent; // Pointer to the meshwork it belongs to.
    size_t _topoIndex; // Index in the meshwork topology.

    unique_ptr<MTriangle> _mTriangle; // pointer to mech triangle

    void updateCoordinate(); // helper function to update coordiante of this triangle

    Compartment* _compartment = nullptr; // The compartment containing this triangle

public:
    mathfunc::Vec< 3, floatingpoint > coordinate; // Coordinate of the center point, updated with updateCoordiante()

    Triangle(Membrane *parent, size_t topoIndex);
    ~Triangle();

    Membrane* getParent()const { return _parent; }
    void setTopoIndex(size_t index) { _topoIndex = index; }
    size_t getTopoIndex() const { return _topoIndex; }

    // Get mech triangle
    MTriangle*       getMTriangle()       { return _mTriangle.get(); }
    const MTriangle* getMTriangle() const { return _mTriangle.get(); }

    /// Get all instances of this class from the SubSystem
    static const vector<Triangle*>& getTriangles() {
        return getElements();
    }

    //@{
    /// SubSystem management, inherited from Trackable
    virtual void addToSubSystem()override { }
    virtual void removeFromSubSystem()override { }
    //@}

    int getType()const;
    void printSelf()const;

    //@{
    /// Implements Movable
    virtual void updatePosition() override;
    //@}
    Compartment* getCompartment() { return _compartment; }


};

#endif
