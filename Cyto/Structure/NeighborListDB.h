//
//  NeighborListDB.h
//  Cyto
//
//  Created by James Komianos on 11/10/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#ifndef __Cyto__NeighborListDB__
#define __Cyto__NeighborListDB__

#include <stdio.h>
#include <list>

#include "common.h"

#include "NeighborList.h"

///NeighborListDB is used to store all NeighborLists in the system

/*! An Object Data Base structure will be used as a container for all main objects: Beads, Filament, Linkers,
 *  Boundary Elements, Motors, and Neighbor Lists. This structure inherits from  list and manage all creations and removing
 *  of objects, as well as some standard list functions and iterators.
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
    
    static NeighborListDB* instance();
    
    ///create a cylinder neighbor list
    CylinderNeighborList* createCylinderNeighborList(float rMax = 0.0, float rMin = 0.0, bool crossFilamentOnly = false) {
        
        CylinderNeighborList* n = new CylinderNeighborList(rMax, rMin, crossFilamentOnly);
        push_back(n);
        
        return n;
    }
    ///create a cylinder neighbor list
    BoundaryElementNeighborList* createBoundaryElementNeighborList(float rMax = 0.0, float rMin = 0.0) {
        
        BoundaryElementNeighborList* n = new BoundaryElementNeighborList(rMax, rMin);
        push_back(n);
        
        return n;
    }

    ///remove a neighborlist
    void removeNeighborList(NeighborList* n) {
        delete n;
        remove(n);
    };
    
    ///reset all neighbors lists
    void resetAll() { for(auto &nlist : *this) nlist->reset(); }
    
    ///add a neighbor to the db. adds to all possible lists
    void addNeighbor(Neighbor* n) { for(auto &nlist : *this) nlist->addNeighbor(n); }
    ///remove a neighbor from the db. removes from all possible lists
    void removeNeighbor(Neighbor* n) { for(auto &nlist : *this) nlist->removeNeighbor(n); }
    
private:
    
    static NeighborListDB* _instance;
    NeighborListDB() {};
};




#endif /* defined(__Cyto__NeighborListDB__) */
