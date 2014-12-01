
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


#ifndef M3SYM_LinkerDB_h
#define M3SYM_LinkerDB_h

#include <list>

#include "common.h"

#include "Linker.h"

//FORWARD DELCARATIONS
class Compartment;
class Cylinder;

/// LinkerDB class is a database for all [Linkers](@ref Linker) in the system
/*!
 *   This LinkerDB inherits from list and manage all creations and removing of
 *   [Linkers](@ref Linker) objects, as well as some standard list functions and iterators.
 *   The [SubSystem] (@ref SubSystem) class calls this database to create and/or remove linkers.
 */
class LinkerDB: private list<Linker*>
{
    typedef list<Linker*> ldb;
    
public:
    using ldb::size;
    using ldb::begin;
    using ldb::end;
    using ldb::erase;
    
    /// Copying is not allowed
    LinkerDB(const LinkerDB &rhs) = delete;
    
    /// Assignment is not allowed
    LinkerDB& operator=(LinkerDB &rhs) = delete;
    
    /// Get instance
    static LinkerDB* instance();
    
    /// Create a new linker
    Linker* createLinker(Cylinder* c1, Cylinder* c2, short linkerType,
                      double position1 = 0.5, double position2 = 0.5, bool creation = false) {
        
        Linker* l = new Linker(c1, c2, linkerType, _currentLinkerID++, position1, position2, creation);
        push_back(l);
        
        return l;
    }

    /// Remove a linker
    void removeLinker(Linker* l) {
        delete l;
        remove(l);
    };
    
private:
    static int _currentLinkerID;   ///< To assign linker IDs
    
    static LinkerDB* _instance;    ///< Singleton instance
    LinkerDB() {};
    
};

#endif
