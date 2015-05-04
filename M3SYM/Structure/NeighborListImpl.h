
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

//FORWARD DECLARATIONS
class Bead;
class Cylinder;
class BoundaryElement;

/// An implementation of NeighborList for Cylinder-Cylinder interactions
class CCNeighborList : public NeighborList {
    
private:
    unordered_map<Cylinder*, vector<Cylinder*>>
        _list; ///< The neighbors list, as a hash map
    
    ///Helper function to update neighbors
    void updateNeighbors(Cylinder* cylinder);
    
public:
    CCNeighborList(float rMax, float rMin=0.0) : NeighborList(rMax, rMin) {}
    
    virtual void addNeighbor(Neighbor* n);
    virtual void removeNeighbor(Neighbor* n);
    
    /// The implementation of this function calls the static version,
    /// all cylinders are dynamic
    virtual void addDynamicNeighbor(Neighbor* n) {addNeighbor(n);}
    virtual void removeDynamicNeighbor(Neighbor* n);
    
    virtual void reset();
    
    /// Get all cylinder neighbors
    vector<Cylinder*> getNeighbors(Cylinder* cylinder);
    
};

/// An implementation of NeighborList for Bead-BoundaryElement interactions
class BBENeighborList : public NeighborList {
    
private:
    unordered_map<BoundaryElement*, vector<Bead*>>
        _list; ///< The neighbors list, as a hash map
    
    ///Helper function to update neighbors
    void updateNeighbors(BoundaryElement* be);
    
public:
    BBENeighborList(float rMax): NeighborList(rMax) {}
    
    virtual void addNeighbor(Neighbor* n);
    virtual void removeNeighbor(Neighbor* n);
    
    virtual void addDynamicNeighbor(Neighbor* n);
    virtual void removeDynamicNeighbor(Neighbor* n);
    
    virtual void reset();
    
    /// Get all Bead neighbors of a boundary element
    vector<Bead*> getNeighbors(BoundaryElement* be);
};

/// An implementation of NeighborList for Cylinder-BoundaryElement interaction
class CBENeighborList : public NeighborList {
    
private:
    unordered_map<BoundaryElement*, vector<Cylinder*>>
        _list; ///< The neighbors list, as a hash map
    
    ///Helper function to update neighbors
    void updateNeighbors(BoundaryElement* be);
    
public:
    CBENeighborList(float rMax): NeighborList(rMax) {}
    
    virtual void addNeighbor(Neighbor* n);
    virtual void removeNeighbor(Neighbor* n);
    
    virtual void addDynamicNeighbor(Neighbor* n);
    virtual void removeDynamicNeighbor(Neighbor* n);
    
    virtual void reset();
    
    /// Get all Cylinder neighbors of a boundary element
    vector<Cylinder*> getNeighbors(BoundaryElement* be);
};

#endif
