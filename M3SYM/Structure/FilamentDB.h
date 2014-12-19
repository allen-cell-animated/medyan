
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

#ifndef M3SYM_FilamentDB_h
#define M3SYM_FilamentDB_h

#include <list>

#include "common.h"

//FORWARD DELCALRATIONS
class Filament;

/// A database for all [Filaments](@ref Filament) in the system.
/*!
 *   This FilamentDB inherits from list and manage all creations and removing of
 *   [Filaments](@ref Filament) objects, as well as some standard list functions 
 *   and iterators.
 */
class FilamentDB: public list<Filament*> {
    typedef list<Filament*> fdb;
    
public:
    using fdb::size;
    using fdb::begin;
    using fdb::end;
    
    /// Copying is not allowed
    FilamentDB(const FilamentDB &rhs) = delete;
    
    /// Assignment is not allowed
    FilamentDB& operator=(FilamentDB &rhs) = delete;
    
    /// Get singleton instance
    static FilamentDB* instance();
    
    /// Add a filament
    void addFilament(Filament* f) { push_back(f); }
    /// Remove a filament
    void removeFilament(Filament* f) { remove(f); };
    
    /// Get current filament ID, and update the ID counter
    int getFilamentID() { return _currentFilamentID++; }
    
private:
    static int _currentFilamentID;  ///< To assign filament ids

    static FilamentDB* _instance;   ///< Singleton instance
    FilamentDB() {};
    
};

#endif
