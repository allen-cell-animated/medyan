//
//  MCylinderDB.h
//  CytoMech
//
//  Created by Konstantin Popov on 7/2/14.
//  Copyright (c) 2014 Konstantin Popov. All rights reserved.
//

#ifndef CytoMech_MCylinderDB_h
#define CytoMech_MCylinderDB_h

#include <iostream>
#include <list>
#include <vector>
#include "Mcommon.h"
#include "MCylinder.h"
#include "MComposite.h"

/*! An Object Data Base structure will be used as a container for all main objects: Beads, Filament, Linkers and Motors. This structure inherits from std:: list and manage all creations and removing of objects, as well as some stabdart list functions and iterators.
 */


class CylinderDB: private std::list<Cylinder*>
{
    typedef std::list<Cylinder*> cdb;
    
public:
    
    //  CylinderDB();
    //  ~CylinderDB();
    
    using cdb::size;
    using cdb::begin;
    using cdb::end;
    using cdb::erase;
    using cdb::remove;
    
    // create new empty cylinder
    Cylinder* CreateCylinder(Filament* pf, Bead* pb) {
        
        Cylinder* c = new Cylinder(pf, pb);
        push_back(c);
        return c ;}
    
    
    
    
    // Remove Cylinder:
    void RemoveCylinder(Cylinder* pc){
        delete pc;
        remove(pc);
        
    }
    
private:
};




#endif
