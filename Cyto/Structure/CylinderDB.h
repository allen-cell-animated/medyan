//
//  CylinderDB.h
//  CytoMech
//
//  Created by Konstantin Popov on 7/2/14.
//  Copyright (c) 2014 Konstantin Popov. All rights reserved.
//

#ifndef CytoMech_CylinderDB_h
#define CytoMech_CylinderDB_h

#include <iostream>
#include <list>
#include <vector>
#include "Cylinder.h"

///Key to access instance of CylinderDB
class CylinderDBKey {friend class Filament; friend class MController; CylinderDBKey(){}; ~CylinderDBKey(){}; };


///CylinderDB class is used to store all Cylinders in the system

/*! An Object Data Base singleton structure will be used as a container for all main objects: Beads, Filament,
 * Linkers, Boundary Elements and Motors. This structure inherits from std:: list and manage all creations
 * and removing of objects, as well as some stabdart list functions and iterators.
 */

class CylinderDB: private std::list<Cylinder*>
{
    typedef std::list<Cylinder*> cdb;
    
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
    
    
    /// get the instance of this singleton
    static CylinderDB* Instance(CylinderDBKey k);
    
    // create new empty cylinder
    Cylinder* CreateCylinder(Filament* pf, Bead* pb, Compartment* c, bool extensionFront = false, bool extensionBack = false) {
        
        Cylinder* pc = new Cylinder(pf, pb, c, extensionFront, extensionBack);
        push_back(pc);
        return pc ;}
    
    
    // Remove Cylinder:
    void RemoveCylinder(Cylinder* pc){
        delete pc;
        remove(pc);
    }
    
private:
    
    static CylinderDB* _instance;
    CylinderDB() {};
};


#endif
