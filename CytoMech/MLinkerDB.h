//
//  MLinkerDB.h
//  CytoMech
//
//  Created by Konstantin Popov on 4/15/14.
//  Copyright (c) 2014 Konstantin Popov. All rights reserved.
//

#ifndef CytoMech_MLinkerDB_h
#define CytoMech_MLinkerDB_h

#include <iostream>
#include <list>
#include "Mcommon.h"
#include "MLinker.h"
#include "MCylinder.h"

/*! An Object Data Base structure will be used as a container for all main objects: Beads, Filament, Linkers and Motors. This structure inherits from std:: list and manage all creations and removing of objects, as well as some stabdart list functions and iterators.
 */

class LinkerDB: private std::list<Linker*>
{
    typedef std::list<Linker*> ldb;
    
public:
    
    
    using ldb::size;
    using ldb::begin;
    using ldb::end;
    using ldb::erase;
    
    
   void CreateLinker(Network* pn, Cylinder* pc1, Cylinder* pc2, double k) {
        
        Linker* pl = new Linker(pn, pc1, pc2, k);
        push_back(pl);
        
        
    }

    
    void RemoveLinker(Linker* pl) {};
};



#endif
