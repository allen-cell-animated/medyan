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

class SpeciesBound;
class Compartment;
class Cylinder;

///Key to access instance of LinkerDB
class LinkerDBKey {friend class SubSystem;
                   friend class LinkerFF;
                   friend class MController;
                   friend class Output;
#ifdef TESTING
                   public:
#endif //TESTING
                   LinkerDBKey() {};
                   ~LinkerDBKey(){}; };

///LinkerDB class is used to store all linkers in the system
/*! An Object Data Base structure will be used as a container for all main objects: Beads, Filament, Linkers 
 *  Boundary Elements, and Motors. This structure inherits from  list and manage all creations and removing
 *  of objects, as well as some stabdart list functions and iterators.
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
    
    static LinkerDB* instance(LinkerDBKey k);
    
    void createLinker(Cylinder* c1, Cylinder* c2, short linkerType, double position1 = 0.5, double position2 = 0.5, bool creation = false) {
        
        Linker* pl = new Linker(c1, c2, linkerType, _currentLinkerID++, position1, position2, creation);
        push_back(pl);
    }

    void removeLinker(Linker* pl) {
        delete pl;
        remove(pl);
    };
    
private:
    ///To assign linker IDs
    static int _currentLinkerID;
    
    static LinkerDB* _instance;
    LinkerDB() {};
    
};

#endif
