//
//  LinkerDB.h
//  CytoMech
//
//  Created by Konstantin Popov on 4/15/14.
//  Copyright (c) 2014 Konstantin Popov. All rights reserved.
//

#ifndef CytoMech_LinkerDB_h
#define CytoMech_LinkerDB_h

#include <iostream>
#include <list>

#include "common.h"
#include "Linker.h"

class Cylinder;

///Key to access instance of LinkerDB
class LinkerDBKey {friend class SubSystem; friend class LinkerFF; LinkerDBKey() {}; ~LinkerDBKey(){}; };

///LinkerDB class is used to store all linkers in the system
/*! An Object Data Base structure will be used as a container for all main objects: Beads, Filament, Linkers 
 *  Boundary Elements, and Motors. This structure inherits from std:: list and manage all creations and removing
 *  of objects, as well as some stabdart list functions and iterators.
 */
class LinkerDB: private std::list<Linker*>
{
    typedef std::list<Linker*> ldb;
    
public:
    using ldb::size;
    using ldb::begin;
    using ldb::end;
    using ldb::erase;
    
    /// Copying is not allowed
    LinkerDB(const LinkerDB &rhs) = delete;
    
    /// Assignment is not allowed
    LinkerDB& operator=(LinkerDB &rhs) = delete;
    
    static LinkerDB* Instance(LinkerDBKey k);
    
    void CreateLinker(Cylinder* pc1, Cylinder* pc2, double k, double position1 = 0.5, double position2 = 0.5) {
        
        Linker* pl = new Linker(pc1, pc2, k, position1, position2);
        push_back(pl);
    }

    void RemoveLinker(Linker* pl) {};
    
private:
    static LinkerDB* _instance;
    LinkerDB() {};
    
};

#endif
