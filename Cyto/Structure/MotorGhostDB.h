//
//  MotorGhostDB.h
//  CytoMech
//
//  Created by Konstantin Popov on 4/16/14.
//  Copyright (c) 2014 Konstantin Popov. All rights reserved.
//

#ifndef CytoMech_MotorGhostDB_h
#define CytoMech_MotorGhostDB_h

#include <iostream>
#include <list>

#include "common.h"
#include "MotorGhost.h"

class Cylinder;

///Key to access instance of MotorGhostDB
class MotorGhostDBKey {friend class SubSystem;
                       friend class MotorGhostFF;
                       friend class MController;
                       friend class Output;
#ifdef TESTING
                       public:
#endif //TESTING
                       MotorGhostDBKey(){}; ~MotorGhostDBKey(){}; };


///MotorGhostDB is used to store all MotorGhosts in the system

/*! An Object Data Base structure will be used as a container for all main objects: Beads, Filament, Linkers, 
 *  Boundary Elements, and Motors. This structure inherits from std:: list and manage all creations and removing 
 *  of objects, as well as some stabdart list functions and iterators.
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
    
    void CreateMotorGhost(Cylinder* pc1, Cylinder* pc2, short motorType, double position1 = 0.5, double position2 = 0.5, bool creation = false) {
    
        MotorGhost* pmg = new MotorGhost(pc1, pc2, motorType, _currentMotorID++, position1, position2, creation);
        push_back(pmg);
        
        //return pmg;
    }
    
    void RemoveMotorGhost(MotorGhost* pmg) {
        delete pmg;
        remove(pmg);
    };
    
private:
    //To assign motor IDs
    static int _currentMotorID;
    
    static MotorGhostDB* _instance;
    MotorGhostDB() {};
};


#endif /* defined(__CytoMech__MotorGhostDB__) */
