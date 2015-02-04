
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

#ifndef M3SYM_BeadDB_h
#define M3SYM_BeadDB_h

#include <list>
#include <vector>

#include "common.h"

//FORWARD DECLARATIONS
class Bead;

/// A database for all [Beads](@ref Bead) in the system.
/*!  
 *   This BeadDB inherits from list and manage all creations and removing of
 *   [Bead](@ref Bead) objects, as well as some standard list functions and iterators.
 */
class BeadDB: private list<Bead*>
{
    typedef list<Bead*> bdb;
    
public:
    using bdb::size;
    using bdb::begin;
    using bdb::end;
    using bdb::erase;
    using bdb::remove;
    
    /// Copying is not allowed
    BeadDB(const BeadDB &rhs) = delete;
    
    /// Assignment is not allowed
    BeadDB& operator=(BeadDB &rhs) = delete;
    
    /// Get the instance of this singleton
    static BeadDB* instance();
    
    /// Add a bead
    void addBead(Bead* b) { push_back(b); }
    /// Remove a bead
    void removeBead(Bead* b) { remove(b); }
    
private:
    static BeadDB* _instance; ///< Singleton instance
    BeadDB() {};
};

#endif 
