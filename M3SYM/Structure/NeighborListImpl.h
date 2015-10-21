
//------------------------------------------------------------------
//  **M3SYM** - Simulation Package for the Mechanochemical
//              Dynamics of Active Networks, 3rd Generation
//
//  Copyright (2014) Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the Papoian lab page for more information:
//  http://papoian.chem.umd.edu/
//------------------------------------------------------------------

#ifndef M3SYM_NeighborListImpl_h
#define M3SYM_NeighborListImpl_h

#include <unordered_map>
#include <vector>

#include "common.h"

#include "NeighborList.h"
#include "DynamicNeighbor.h"

//FORWARD DECLARATIONS
class Cylinder;
class Bubble;
class BoundaryElement;

/// An implementation of NeighborList for Cylinder-Cylinder interactions
/// This can be a half or full list depending on the usage.
class CylinderCylinderNL : public NeighborList {
    
private:
    unordered_map<Cylinder*, vector<Cylinder*>> _list;
    ///< The neighbors list, as a hash map
    
    bool _full; ///<Specifying whether this is a full or half list
    
    ///Helper function to update neighbors
    ///@param runtime - specifying whether the cylinder is being
    ///created/destroyed at runtime vs at a full neighbor list update.
    void updateNeighbors(Cylinder* cylinder, bool runtime = false);
    
public:
    CylinderCylinderNL(float rMax, float rMin=0.0, bool full = false)
    
        : NeighborList(rMax, rMin), _full(full) {}
    
    virtual void addNeighbor(Neighbor* n);
    virtual void removeNeighbor(Neighbor* n);
    
    //@{
    /// The implementation of these functions calls the static version,
    /// all cylinders are dynamic
    virtual void addDynamicNeighbor(DynamicNeighbor* n) {addNeighbor(n);}
    virtual void removeDynamicNeighbor(DynamicNeighbor* n) {removeNeighbor(n);}
    //@}
    
    virtual void reset();
    
    /// Get all cylinder neighbors
    vector<Cylinder*> getNeighbors(Cylinder* cylinder);
    
};

/// An implementation of NeighborList for BoundaryElement-Cylinder interactions
class BoundaryCylinderNL : public NeighborList {
    
private:
    unordered_map<BoundaryElement*, vector<Cylinder*>> _list;
    ///< The neighbors list, as a hash map
    
    ///Helper function to update neighbors
    void updateNeighbors(BoundaryElement* be);
    
public:
    BoundaryCylinderNL(float rMax): NeighborList(rMax) {}
    
    virtual void addNeighbor(Neighbor* n);
    virtual void removeNeighbor(Neighbor* n);
    
    virtual void addDynamicNeighbor(DynamicNeighbor* n);
    virtual void removeDynamicNeighbor(DynamicNeighbor* n);
    
    virtual void reset();
    
    /// Get all Cylinder neighbors of a boundary element
    vector<Cylinder*> getNeighbors(BoundaryElement* be);
};


/// An implementation of NeighborList for BoundaryElement-Bubble interactions
class BoundaryBubbleNL : public NeighborList {
    
private:
    unordered_map<BoundaryElement*, vector<Bubble*>> _list;
    ///< The neighbors list, as a hash map
    
    ///Helper function to update neighbors
    void updateNeighbors(BoundaryElement* be);
    
public:
    BoundaryBubbleNL(float rMax): NeighborList(rMax) {}
    
    virtual void addNeighbor(Neighbor* n);
    virtual void removeNeighbor(Neighbor* n);
    
    virtual void addDynamicNeighbor(DynamicNeighbor* n);
    virtual void removeDynamicNeighbor(DynamicNeighbor* n);
    
    virtual void reset();
    
    /// Get all Bubble neighbors of a boundary element
    vector<Bubble*> getNeighbors(BoundaryElement* be);
};

/// An implementation of NeighborList for Bubble-Bubble interactions
/// @note - This is currently implemented as a half list only
class BubbleBubbleNL : public NeighborList {
    
private:
    unordered_map<Bubble*, vector<Bubble*>> _list;
    ///< The neighbors list, as a hash map
    
    ///Helper function to update neighbors
    void updateNeighbors(Bubble* bb);
    
public:
    BubbleBubbleNL(float rMax): NeighborList(rMax) {}
    
    virtual void addNeighbor(Neighbor* n);
    virtual void removeNeighbor(Neighbor* n);
    
    //@{
    /// The implementation of these functions calls the static version,
    /// all Bubbles are dynamic
    virtual void addDynamicNeighbor(DynamicNeighbor* n) {addNeighbor(n);}
    virtual void removeDynamicNeighbor(DynamicNeighbor* n) {removeNeighbor(n);}
    //@}
    
    virtual void reset();
    
    /// Get all Bubble neighbors of a bubble
    vector<Bubble*> getNeighbors(Bubble* bb);
};

/// An implementation of NeighborList for Bubble-Cylinder interactions
class BubbleCylinderNL : public NeighborList {
    
private:
    unordered_map<Bubble*, vector<Cylinder*>> _list;
    ///< The neighbors list, as a hash map
    
    ///Helper function to update neighbors
    void updateNeighbors(Bubble* bb);
    
public:
    BubbleCylinderNL(float rMax): NeighborList(rMax) {}
    
    virtual void addNeighbor(Neighbor* n);
    virtual void removeNeighbor(Neighbor* n);
    
    //@{
    /// The implementation of these functions calls the static version,
    /// all Bubbles and Cylinders are dynamic
    virtual void addDynamicNeighbor(DynamicNeighbor* n) {addNeighbor(n);}
    virtual void removeDynamicNeighbor(DynamicNeighbor* n) {removeNeighbor(n);}
    //@}
    
    virtual void reset();
    
    /// Get all Cylinder neighbors of a bubble
    vector<Cylinder*> getNeighbors(Bubble* bb);
};





#endif
