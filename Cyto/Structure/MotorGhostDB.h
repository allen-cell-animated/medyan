
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

#include "MotorGhost.h"

//FORWARD DECLARATIONS
class Cylinder;

/// MotorGhostDB class is a database for all [MotorGhosts](@ref MotorGhost) in the system
/*!
 *   This MotorGhostDB inherits from list and manage all creations and removing of
 *   [MotorGhosts](@ref MotorGhost) objects, as well as some standard list functions and iterators.
 *   The [SubSystem] (@ref SubSystem) class calls this database to create and/or remove MotorGhosts.
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
    
    /// Create a motor ghost
    MotorGhost* createMotorGhost(Cylinder* c1, Cylinder* c2, short motorType,
                                 double position1 = 0.5, double position2 = 0.5, bool creation = false) {
    
        MotorGhost* mg = new MotorGhost(c1, c2, motorType, _currentMotorID++, position1, position2, creation);
        push_back(mg);
        
        return mg;
    }
    
    /// Remove a motor ghost
    void removeMotorGhost(MotorGhost* mg) {
        delete mg;
        remove(mg);
    };
    
private:
    static int _currentMotorID; ///< To assign motor IDs
    
    static MotorGhostDB* _instance; ///< Singleton instance
    MotorGhostDB() {};
};


#endif
