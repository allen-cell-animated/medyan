//
//  MMotorGhostDB.h
//  CytoMech
//
//  Created by Konstantin Popov on 4/16/14.
//  Copyright (c) 2014 Konstantin Popov. All rights reserved.
//

#ifndef CytoMech_MMotorGhostDB_h
#define CytoMech_MMotorGhostDB_h

#include <iostream>
#include <list>
#include "Mcommon.h"
#include "MMotorGhost.h"
#include "MCylinder.h"

///Key to access instance of MotorGhostDB
class MotorGhostDBKey {friend class System; MotorGhostDBKey(); ~MotorGhostDBKey(); };

/*! An Object Data Base structure will be used as a container for all main objects: Beads, Filament, Linkers and Motors. This structure inherits from std:: list and manage all creations and removing of objects, as well as some stabdart list functions and iterators.
 */

class MotorGhostDB: private std::list<MotorGhost*>
{
    typedef std::list<MotorGhost*> mgdb;
    
public:
    
    
    using mgdb::size;
    using mgdb::begin;
    using mgdb::end;
    using mgdb::erase;
    
    /// Copying is not allowed
    MotorGhostDB(const MotorGhostDB &rhs) = delete;
    
    /// Assignment is not allowed
    MotorGhostDB& operator=(MotorGhostDB &rhs) = delete;
    
    static MotorGhostDB* Instance(MotorGhostDBKey k);
    
    void CreateMotorGhost(Network* pn, Cylinder* pc1, Cylinder* pc2, double k, double position1, double position2) {
    
        MotorGhost* pmg = new MotorGhost(pn, pc1, pc2, k, position1, position2);
        push_back(pmg);
        
        //return pmg;
    }
    
    
    void RemoveMotorGhost(MotorGhost* pmg) {};
    
private:
    static MotorGhostDB* _instance;
    MotorGhostDB() {};
};


#endif
