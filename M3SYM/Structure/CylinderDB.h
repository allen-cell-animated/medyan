
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

#ifndef M3SYM_CylinderDB_h
#define M3SYM_CylinderDB_h

#include <list>
#include <vector>

#include "common.h"

#include "Cylinder.h"

/// CylinderDB class is a database for all [Cylinders](@ref Cylinder) in the system
/*!
 *   This CylinderDB inherits from list and manage all creations and removing of
 *   [Cylinders](@ref Cylinder) objects, as well as some standard list functions and iterators.
 *   The [Filament] (@ref Filament) class calls this database to create and/or remove cylinders.
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
    
    /// Create new empty cylinder
    Cylinder* createCylinder(Filament* pf, Bead* firstBead, Bead* secondBead, int positionFilament, 
                             bool extensionFront = false, bool extensionBack = false, bool creation = false) {
        
        Cylinder* pc = new Cylinder(pf, firstBead, secondBead, positionFilament, _currentCylinderID++,
                                    extensionFront, extensionBack, creation);
        push_back(pc);
        return pc;
    }
    
    
    /// Remove Cylinder
    void removeCylinder(Cylinder* c){
        delete c;
        remove(c);
    }
    
private:
    int _currentCylinderID; ///< Current running cylinder ID
    
    static CylinderDB* _instance; ///< Singleton instance
    CylinderDB() {};
};


#endif
