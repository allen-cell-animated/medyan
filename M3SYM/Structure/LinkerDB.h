
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

//FORWARD DELCARATIONS
class Linker;

/// LinkerDB class is a database for all [Linkers](@ref Linker) in the system
/*!
 *   This LinkerDB inherits from list and manage all creations and removing of
 *   [Linkers](@ref Linker) objects, as well as some standard list functions and iterators.
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
    void addLinker(Linker* l) { push_back(l); }
    /// Remove a linker
    void removeLinker(Linker* l) { remove(l); };
    
    /// Get current linker ID, and update the ID counter
    int getLinkerID() { return _currentLinkerID++; }
    
private:
    static int _currentLinkerID;   ///< To assign linker IDs
    
    static LinkerDB* _instance;    ///< Singleton instance
    LinkerDB() {};
    
};

#endif
