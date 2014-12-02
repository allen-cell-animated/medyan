
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

#ifndef M3SYM_MotorGhostDB_h
#define M3SYM_MotorGhostDB_h

#include <list>

#include "common.h"

//FORWARD DECLARATIONS
class MotorGhost;

/// MotorGhostDB class is a database for all [MotorGhosts](@ref MotorGhost) in the system
/*!
 *   This MotorGhostDB inherits from list and manage all creations and removing of
 *   [MotorGhosts](@ref MotorGhost) objects, as well as some standard list functions and iterators.
 */
class MotorGhostDB: private list<MotorGhost*>
{
    typedef list<MotorGhost*> mgdb;
    
public:
    using mgdb::size;
    using mgdb::begin;
    using mgdb::end;
    using mgdb::erase;
    
    /// Copying is not allowed
    MotorGhostDB(const MotorGhostDB &rhs) = delete;
    
    /// Assignment is not allowed
    MotorGhostDB& operator=(MotorGhostDB &rhs) = delete;
    
    /// Get instance
    static MotorGhostDB* instance();
    
    /// Add a motor ghost
    void addMotorGhost(MotorGhost* mg) {push_back(mg); }
    /// Remove a motor ghost
    void removeMotorGhost(MotorGhost* mg) { remove(mg); }
    
    /// Get current motor ID, and update the ID counter
    int getMotorID() { return _currentMotorID++; }
    
private:
    static int _currentMotorID; ///< To assign motor IDs
    
    static MotorGhostDB* _instance; ///< Singleton instance
    MotorGhostDB() {};
};


#endif
