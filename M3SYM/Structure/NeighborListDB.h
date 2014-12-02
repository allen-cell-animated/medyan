
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

#ifndef M3SYM_NeighborListDB_h
#define M3SYM_NeighborListDB_h

#include <list>

#include "common.h"

#include "NeighborList.h"

//FORWARD DELCARATIONS
class Neighbor;

/// A database to store all [NeighborLists](@ref NeighborList) in the system.
/*!
 *   This NeighborListDB inherits from list and manage all creations and removing of
 *   [NeighborLists](@ref NeighborList) objects, as well as some standard list functions and iterators.
 */
class NeighborListDB: private list<NeighborList*>
{
    typedef list<NeighborList*> ndb;
    
public:
    using ndb::size;
    using ndb::begin;
    using ndb::end;
    using ndb::erase;
    
    /// Copying is not allowed
    NeighborListDB(const NeighborListDB &rhs) = delete;
    
    /// Assignment is not allowed
    NeighborListDB& operator=(NeighborListDB &rhs) = delete;
    
    /// Get instance
    static NeighborListDB* instance();
    
    /// Add a cylinder neighbor list
    void addNeighborList(NeighborList* n) { push_back(n); }

    /// Remove a neighborlist
    void removeNeighborList(NeighborList* n) { remove(n); };
    
    /// Reset all neighbors lists
    void resetAll() { for(auto &nlist : *this) nlist->reset(); }
    
    /// Add a neighbor to the db. adds to all possible lists
    void addNeighbor(Neighbor* n) { for(auto &nlist : *this) nlist->addNeighbor(n); }
    /// Remove a neighbor from the db. removes from all possible lists
    void removeNeighbor(Neighbor* n) { for(auto &nlist : *this) nlist->removeNeighbor(n); }
    
private:
    static NeighborListDB* _instance;   ///< Singleton instance
    NeighborListDB() {};
};

#endif
