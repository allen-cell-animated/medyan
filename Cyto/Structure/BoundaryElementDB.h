//
//  BoundaryElementDB.h
//  Cyto
//
//  Created by James Komianos on 9/2/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#ifndef __Cyto__BoundaryElementDB__
#define __Cyto__BoundaryElementDB__

#include "common.h"
#include "BoundaryElementImpl.h"
#include <iostream>

///Key to access instance of BoundaryElementDB
class BoundaryElementDBKey {friend class BoundarySurface; friend class Bead; friend class BoundaryFF; BoundaryElementDBKey(){}; public: ~BoundaryElementDBKey(){};};


/*! An Object Data Base singleton structure will be used as a container for all main objects: Boundary Elements, Beads, Filament, Linkers and Motors. This structure inherits from std:: list and manage all creations and removing of objects, as well as some stabdart list functions and iterators.
 */


class BoundaryElementDB: private std::list<BoundaryElement*>
{
    typedef std::list<BoundaryElement*> bedb;
    
public:
    using bedb::size;
    using bedb::begin;
    using bedb::end;
    using bedb::erase;
    using bedb::remove;
    
    /// Copying is not allowed
    BoundaryElementDB(const BoundaryElementDB &rhs) = delete;
    
    /// Assignment is not allowed
    BoundaryElementDB& operator=(BoundaryElementDB &rhs) = delete;
    
    /// get the instance of this singleton
    static BoundaryElementDB* Instance(BoundaryElementDBKey k);
    
//    BoundaryElement* CreateTriangleBoundaryElement(std::vector<double> coords, float size) {
//        
//        //BoundaryElement* b = new TriangleBoundaryElement(coords, size);
//        //push_back(b);
//        //return b ;
//        return nullptr;
//    }
    
    /// create a new plane boundary element
    BoundaryElement* CreatePlaneBoundaryElement(std::vector<double>& coords, std::vector<double>& normal) {
        
        BoundaryElement* b = new PlaneBoundaryElement(coords, normal);
        push_back(b);
        return b ;
    }
    
    // Remove boundary element
    void RemoveBoundaryElement(BoundaryElement* pb){
        
        remove(pb);
        delete pb;
        
    }
private:
    static BoundaryElementDB* _instance;
    BoundaryElementDB() {};
    
};








#endif /* defined(__Cyto__BoundaryElementDB__) */
