
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

#ifndef M3SYM_CylinderDB_h
#define M3SYM_CylinderDB_h

#include <list>
#include <vector>

#include "common.h"

//FORWARD DECLARATIONS
class Cylinder;

/// A database for all [Cylinders](@ref Cylinder) in the system.
/*!
 *   This CylinderDB inherits from list and manage all creations and removing of
 *   [Cylinders](@ref Cylinder) objects, as well as some standard list functions and 
 *   iterators.
 */
class CylinderDB: private list<Cylinder*>
{
    typedef list<Cylinder*> cdb;
    
public:
    using cdb::size;
    using cdb::begin;
    using cdb::end;
    using cdb::erase;
    using cdb::remove;
    
    /// Copying is not allowed
    CylinderDB(const CylinderDB &rhs) = delete;
    
    /// Assignment is not allowed
    CylinderDB& operator=(CylinderDB &rhs) = delete;
    
    /// Get the instance of this singleton
    static CylinderDB* instance();
    
    /// Add a cylinder
    void addCylinder(Cylinder* c){ push_back(c); }
    /// Remove cylinder
    void removeCylinder(Cylinder* c){ remove(c); }
    
    /// Get current cylinder ID, and update the ID counter
    int getCylinderID() { return _currentCylinderID++; }
    
private:
    static int _currentCylinderID; ///< Current running cylinder ID
    
    static CylinderDB* _instance; ///< Singleton instance
    CylinderDB() {};
};


#endif
