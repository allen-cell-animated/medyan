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

#include "common.h"
#include "Cylinder.h"

///Key to access instance of CylinderDB
class CylinderDBKey {friend class Filament;
                     friend class Controller;
                     friend class MController;
                     friend class VolumeCylindricalFF;
#ifdef TESTING
                     public:
#endif //TESTING
                     CylinderDBKey(){}; ~CylinderDBKey(){}; };


///CylinderDB class is used to store all Cylinders in the system

/*! An Object Data Base singleton structure will be used as a container for all main objects: Beads, Filament,
 * Linkers, Boundary Elements, Motors, and Neighbor Lists. This structure inherits from  list and manage all creations
 * and removing of objects, as well as some standard list functions and iterators.
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
    
    
    /// get the instance of this singleton
    static CylinderDB* instance(CylinderDBKey k);
    
    // create new empty cylinder
    Cylinder* createCylinder(Filament* pf, Bead* firstBead, Bead* secondBead, int positionFilament, 
                             bool extensionFront = false, bool extensionBack = false, bool creation = false) {
        
        Cylinder* pc = new Cylinder(pf, firstBead, secondBead, positionFilament, _currentCylinderID++,
                                    extensionFront, extensionBack, creation);
        push_back(pc);
        return pc;
    }
    
    
    // Remove Cylinder:
    void removeCylinder(Cylinder* c){
        delete c;
        remove(c);
    }
    
private:
    int _currentCylinderID; ///current running cylinder ID
    
    static CylinderDB* _instance;
    CylinderDB() {};
};


#endif
