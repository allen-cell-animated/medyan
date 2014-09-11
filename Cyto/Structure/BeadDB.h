//
//  BeadDB.h
//  CytoMech
//
//  Created by Konstantin Popov on 4/15/14.
//  Copyright (c) 2014 Konstantin Popov. All rights reserved.
//

#ifndef CytoMech_BeadDB_h
#define CytoMech_BeadDB_h

#include <iostream>
#include <list>
#include <vector>
#include "Bead.h"

///Key to access instance of BeadDB
class BeadDBKey {friend class Cylinder;
                 friend class Filament;
                 friend class CGMethod;
                 BeadDBKey() {};
                 public: ~BeadDBKey() {}; };


///BeadDB class is a database for all beads in the system

/*!  An Object Data Base singleton structure will be used as a container for all main objects: Beads, Filament, Linkers
 *   Boundary Elements and Motors. This structure inherits from std:: list and manage all creations and removing of 
 *   objects, as well as some stabdart list functions and iterators.
 */
class BeadDB: private std::list<Bead*>
{
    typedef std::list<Bead*> bdb;
    
public:
    
    using bdb::size;
    using bdb::begin;
    using bdb::end;
    using bdb::erase;
    using bdb::remove;
    
    /// Copying is not allowed
    BeadDB(const BeadDB &rhs) = delete;
    
    /// Assignment is not allowed
    BeadDB& operator=(BeadDB &rhs) = delete;
    
    /// get the instance of this singleton
    static BeadDB* Instance(BeadDBKey k);
    
    /// create a new bead with no ccordinates and
    Bead* CreateBead() {
        
        Bead* b = new Bead();
        push_back(b);
        return b ;}
    
    /// Create bead with a given coordinate on a given filament:
    Bead* CreateBead(std::vector<double> v) {
        
        Bead* b = new Bead(v);
        
        push_back(b);
        std::cout<<"bead created"<<std::endl;
        return b ;}
    
    
    /// Remove bead:
    void RemoveBead(Bead* pb){
        ///Also need to clean all neighbours lists!
        delete pb;
        remove(pb);
        
    }
private:
    static BeadDB* _instance;
    BeadDB() {};
};

#endif /* defined(__CytoMech__BeadDB__) */
