
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

#ifndef M3SYM_BranchingPointDB_h
#define M3SYM_BranchingPointDB_h

#include <list>

#include "common.h"

//FORWARD DECLARATIONS
class BranchingPoint;

/// A database for all [BranchingPoints](@ref BranchingPoint) in the system.
/*!
 *   This BranchingPointDB inherits from list and manage all creations and removing of
 *   [BranchingPoints](@ref BranchingPoint) objects, as well as some standard list functions and iterators.
 */
class BranchingPointDB: private list<BranchingPoint*>
{
    typedef list<BranchingPoint*> bpdb;
    
public:
    using bpdb::size;
    using bpdb::begin;
    using bpdb::end;
    using bpdb::erase;
    
    /// Copying is not allowed
    BranchingPointDB(const BranchingPointDB &rhs) = delete;
    
    /// Assignment is not allowed
    BranchingPointDB& operator=(BranchingPointDB &rhs) = delete;
    
    /// Get instance
    static BranchingPointDB* instance();
    
    /// Create a new branch point
    void addBranchingPoint(BranchingPoint* b) { push_back(b); }
    /// Remove a branch point
    void removeBranchingPoint(BranchingPoint* b) { remove(b); };
    
    /// Get current branch point ID, and update the ID counter
    int getBranchID() { return _currentBranchID++; }
    
private:
    static int _currentBranchID;   ///< To assign branch IDs
    
    static BranchingPointDB* _instance;    ///< Singleton instance
    BranchingPointDB() {};
};


#endif
