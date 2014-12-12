
//------------------------------------------------------------------
//  **M3SYM** - Simulation Package for the Mechanochemical
//              Dynamics of Active Networks, 3rd Generation
//
//  Copyright (2014) Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the Papoian lab page for installation and documentation:
//  http://papoian.chem.umd.edu/
//------------------------------------------------------------------

#ifndef M3SYM_NeighborList_h
#define M3SYM_NeighborList_h

#include <unordered_map>
#include <vector>

#include "common.h"

//FORWARD DECLARATIONS
class Neighbor;
class Cylinder;
class BoundaryElement;
class Bead;

/// To hold an external neighbor list of general type.

/*!
 *  This class is used to hold any neighbor list. Contains a map of neighbors as well as
 *  min and max cutoffs for generation of the list. This class is abstract and must be
 *  implemented by writing functionality to add and update a neighbor.
 * 
 *  The neighbor list contains a function to reset, which uses the databases to clear
 *  and update the list.
 */
class NeighborList {
    
protected:
    unordered_map<Neighbor*, vector<Neighbor*>> _list; ///< The neighbors list, as a hash map
    float _rMax;  ///< max distance cutoff
    float _rMin;  ///< min distance cutoff
    
public:
    ///Constructor and destructor
    NeighborList(float rMax = 0.0, float rMin = 0.0) : _rMax(rMax), _rMin(rMin) {}
    
    ///Destructor
    /// @note noexcept is important here. Otherwise, gcc flags the constructor as potentially throwing,
    /// which in turn disables move operations by the STL containers. This behaviour is a gcc bug
    /// (as of gcc 4.703), and will presumbaly be fixed in the future.
    virtual ~NeighborList() noexcept {}
    
    /// Remove a neighbor if possible
    virtual void removeNeighbor(Neighbor* n) {_list.erase(n);}
    
    /// Re-initialize the neighborlist
    virtual void reset() {
        //loop through all neighbor keys
        for(auto it = _list.begin(); it != _list.end(); it++) {
            
            it->second.clear(); ///clear vector of neighbors
            updateNeighbors(it->first);
        }
    }
    
    /// Add neighbor
    virtual void addNeighbor(Neighbor* n) = 0;
    /// Update a neighbor
    virtual void updateNeighbors(Neighbor* n) = 0;
    
    /// Add a dynamic neighbor to the system.
    /// For BoundaryElementNeighborList list, this will be Bead.
    /// For CylinderNeighborList, all elements in neighbors list are dynamic.
    virtual void addDynamicNeighbor(Neighbor* n) = 0;
    
    /// Remove a dynamic neighbor from the system.
    /// For BoundaryElementNeighborList list, this will be Bead.
    /// For CylinderNeighborList, all elements in neighbors list are dynamic.
    virtual void removeDynamicNeighbor(Neighbor* n) = 0;
};


/// An implementation of NeighborList for Cylinder-Cylinder interactions
class CylinderNeighborList : public NeighborList {
    
private:
    bool _crossFilamentOnly; ///< Whether to include cylinders in same filament
    
public:
    CylinderNeighborList(float rMax = 0.0, float rMin = 0.0, bool crossFilamentOnly = false)
                         : NeighborList(rMax, rMin), _crossFilamentOnly(crossFilamentOnly) {}
    
    virtual void addNeighbor(Neighbor* n);
    virtual void updateNeighbors(Neighbor* n);
    
    /// The implementation of this function calls the static version, all cylinders are dynamic
    virtual void addDynamicNeighbor(Neighbor* n) {addNeighbor(n);}
    virtual void removeDynamicNeighbor(Neighbor* n);
    
    /// Get all cylinder neighbors
    vector<Cylinder*> getNeighbors(Cylinder* cylinder);

};

/// An implementation of NeighborList for Bead-BoundaryElement interactions
class BoundaryElementNeighborList : public NeighborList {
    
public:
    BoundaryElementNeighborList(float rMax, float rMin) : NeighborList(rMax, rMin) {}

    virtual void addNeighbor(Neighbor* n);
    virtual void updateNeighbors(Neighbor* n);
    
    virtual void addDynamicNeighbor(Neighbor* n);
    virtual void removeDynamicNeighbor(Neighbor* n);
    
    /// Get all Bead neighbors of a boundary element
    vector<Bead*> getNeighbors(BoundaryElement* be);
};



#endif
