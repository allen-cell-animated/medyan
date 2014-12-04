
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

#ifndef M3SYM_BranchPointDB_h
#define M3SYM_BranchPointDB_h

#include <list>

#include "common.h"

//FORWARD DECLARATIONS
class BranchPoint;

/// A database for all [BranchPoints](@ref BranchPoint) in the system.
/*!
 *   This BranchPointDB inherits from list and manage all creations and removing of
 *   [BranchPoints](@ref BranchPoint) objects, as well as some standard list functions and iterators.
 */
class BranchPointDB: private list<BranchPoint*>
{
    typedef list<BranchPoint*> bpdb;
    
public:
    using bpdb::size;
    using bpdb::begin;
    using bpdb::end;
    using bpdb::erase;
    
    /// Copying is not allowed
    BranchPointDB(const BranchPointDB &rhs) = delete;
    
    /// Assignment is not allowed
    BranchPointDB& operator=(BranchPointDB &rhs) = delete;
    
    /// Get instance
    static BranchPointDB* instance();
    
    /// Create a new branch point
    void addBranchPoint(BranchPoint* b) { push_back(b); }
    /// Remove a branch point
    void removeBranchPoint(BranchPoint* b) { remove(b); };
    
    /// Get current branch point ID, and update the ID counter
    int getBranchID() { return _currentBranchID++; }
    
private:
    static int _currentBranchID;   ///< To assign branch IDs
    
    static BranchPointDB* _instance;    ///< Singleton instance
    BranchPointDB() {};
};


#endif
