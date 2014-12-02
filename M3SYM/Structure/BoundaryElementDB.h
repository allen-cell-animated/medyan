
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

#ifndef M3SYM_BoundaryElementDB_h
#define M3SYM_BoundaryElementDB_h

#include <list>

#include "common.h"

//FORWARD DECLARATIONS
class BoundaryElement;

/// A database for all [BoundaryElements](@ref BoundaryElement) in the system.
/*!
 *   This BoundaryElementDB inherits from list and manage all creations and removing of
 *   [BoundaryElement](@ref BoundaryElement) objects, as well as some standard list functions and iterators.
 */
class BoundaryElementDB: private list<BoundaryElement*>
{
    typedef list<BoundaryElement*> bedb;
    
public:
    using bedb::size;
    using bedb::begin;
    using bedb::end;
    using bedb::erase;
    using bedb::remove;
    
    /// Copying is not allowed
    BoundaryElementDB(const BoundaryElementDB &rhs) = delete;
    
    /// Assignment is not allowed
    BoundaryElementDB& operator=(BoundaryElementDB &rhs) = delete;
    
    /// Get the instance of this singleton
    static BoundaryElementDB* instance();
    
    /// Add a boundary element
    void addBoundaryElement(BoundaryElement* b) { push_back(b); }
    /// Remove boundary element
    void removeBoundaryElement(BoundaryElement* b){ remove(b); }
    
private:
    static BoundaryElementDB* _instance; ///< Singleton instance
    BoundaryElementDB() {};
    
};


#endif
