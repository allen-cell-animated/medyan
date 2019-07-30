#ifndef MEDYAN_Structure_SurfaceMesh_Edge_Hpp
#define MEDYAN_Structure_SurfaceMesh_Edge_Hpp

/*
 
An unordered edge contains 2 vertices.
 
*/

#include "common.h"
#include "Database.h"
#include "Trackable.h"
#include "DynamicNeighbor.h"
#include "MathFunctions.h"
#include "Movable.h"
#include "Structure/CellList.hpp"

// Forward declarations
class Compartment;
class Membrane;

class Edge:
    public Movable,
    public Trackable,
    public Database< Edge, false > {

private:

    Membrane* _parent;
    size_t _topoIndex; // Index in the meshwork topology.

    void updateCoordinate(); // helper function to update coordiante of this edge

    cell_list::CellListElementUser< Edge, Compartment > _cellElement;

public:
    Edge(Membrane *parent, size_t topoIndex);
    ~Edge();

    Membrane* getParent()const { return _parent; }
    void setTopoIndex(size_t index) { _topoIndex = index; }

    mathfunc::Vec< 3, floatingpoint > coordinate; // Coordinate of the mid point, updated with updateCoordiante()

    /// Get all instances of this class from the SubSystem
    static const vector<Edge*>& getEdges() {
        return getElements();
    }

    int getType()const;
    void printSelf()const;

    //@{
    /// SubSystem management, inherited from Trackable
    virtual void addToSubSystem()override { }
    virtual void removeFromSubSystem()override { }
    //@}

    //@{
    /// Implements Movable
    virtual void updatePosition() override;
    //@}
    Compartment* getCompartment() const { return _cellElement.manager->getHeadPtr(_cellElement); }

};

#endif
