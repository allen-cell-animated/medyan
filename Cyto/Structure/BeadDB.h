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

#include "common.h"
#include "Bead.h"

///BeadDB class is a database for all beads in the system

/*!  An Object Data Base singleton structure will be used as a container for all main objects: Beads, Filament, Linkers
 *   Boundary Elements and Motors. This structure inherits from list and manage all creations and removing of
 *   objects, as well as some standard list functions and iterators.
 */
class BeadDB: private list<Bead*>
{
    typedef list<Bead*> bdb;
    
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
    static BeadDB* instance();
    
    /// create a new bead with no ccordinates and
    Bead* createBead(int positionFilament) {
        
        Bead* b = new Bead(positionFilament);
        push_back(b);
        return b ;}
    
    /// Create bead with a given coordinate on a given filament:
    Bead* createBead(vector<double>& v, int positionFilament) {
        
        Bead* b = new Bead(v, positionFilament);
        push_back(b);
        return b ;}
    
    
    /// Remove bead:
    void removeBead(Bead* b){
        delete b;
        remove(b);
        
    }
private:
    static BeadDB* _instance;
    BeadDB() {};
};

#endif /* defined(__CytoMech__BeadDB__) */
