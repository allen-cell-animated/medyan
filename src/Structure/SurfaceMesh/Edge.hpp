#ifndef MEDYAN_Structure_SurfaceMesh_Edge_Hpp
#define MEDYAN_Structure_SurfaceMesh_Edge_Hpp

/*
 
An unordered edge contains 2 vertices.
 
*/

#include "common.h"
#include "MathFunctions.h"
#include "Structure/CellList.hpp"
#include "Structure/Database.h"
#include "Structure/DynamicNeighbor.h"
#include "Structure/Movable.h"
#include "Structure/Trackable.h"

// Forward declarations
class Compartment;
class Membrane;

class Edge:
    public Movable,
    public Trackable,
    public Database< Edge, false > {

private:

    Membrane* parent_;
    size_t topoIndex_; // Index in the meshwork topology.

    void updateCoordinate(); // helper function to update coordiante of this edge

    cell_list::CellListElementUser< Edge, Compartment > _cellElement;

public:
    Edge(Membrane *parent, size_t topoIndex);
    ~Edge();

    Membrane* getParent()const { return parent_; }
    void setTopoIndex(size_t index) { topoIndex_ = index; }

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
