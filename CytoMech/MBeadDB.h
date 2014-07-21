//
//  MBeadDB.h
//  CytoMech
//
//  Created by Konstantin Popov on 4/15/14.
//  Copyright (c) 2014 Konstantin Popov. All rights reserved.
//

#ifndef CytoMech_MBeadDB_h
#define CytoMech_MBeadDB_h

#include <iostream>
#include <list>
#include <vector>
#include "MBead.h"
#include "MComposite.h"

/*! An Object Data Base structure will be used as a container for all main objects: Beads, Filament, Linkers and Motors. This structure inherits from std:: list and manage all creations and removing of objects, as well as some stabdart list functions and iterators.
 */


class BeadDB: private std::list<Bead*>
{
    typedef std::list<Bead*> bdb;
    
public:
    
    //  BeadDB();
    //  ~BeadDB();
    
    using bdb::size;
    using bdb::begin;
    using bdb::end;
    using bdb::erase;
    using bdb::remove;
    
    // create a new bead with no ccordinates and
    Bead* CreateBead() {
        
        Bead* b = new Bead();
        push_back(b);
        return b ;}
    
    // Create  bid with a given coordinate on a given filament:
    Bead* CreateBead(std::vector<double> v) {
        
        
        Bead* b = new Bead(v);
        
        push_back(b);
        std::cout<<"bead created"<<std::endl;
        return b ;}
    
    
    // Remove bead:
    void RemoveBead(Bead* pb){
        pb->getParent()->DeleteBead(pb);
        ///Also need to clean all neighbours lists!
        delete pb;
        remove(pb);
        
    }
    
private:
};





#endif
