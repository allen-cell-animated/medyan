//
//  BoundaryElementDB.h
//  Cyto
//
//  Created by James Komianos on 9/2/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#ifndef __Cyto__BoundaryElementDB__
#define __Cyto__BoundaryElementDB__

#include <iostream>

#include "common.h"
#include "BoundaryElementImpl.h"

///Key to access instance of BoundaryElementDB
class BoundaryElementDBKey {friend class BoundarySurface;
                            friend class Bead;
                            friend class BoundaryFF;
#ifdef TESTING
                            public:
#endif //TESTING
                            BoundaryElementDBKey(){}; public: ~BoundaryElementDBKey(){};};


/*! An Object Data Base singleton structure will be used as a container for all main objects: 
 *  Boundary Elements, Beads, Filament, Linkers and Motors. This structure inherits from 
 *  list and manage all creations and removing of objects, as well as some standard
 *  list functions and iterators.
 */


class BoundaryElementDB: private list<BoundaryElement*>
{
    typedef list<BoundaryElement*> bedb;
    
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
    static BoundaryElementDB* instance(BoundaryElementDBKey k);
    
    /// create a new plane boundary element
    BoundaryElement* createPlaneBoundaryElement(vector<double>& coords, vector<double>& normal, double repulsConst, double screenLength) {
        
        BoundaryElement* b = new PlaneBoundaryElement(coords, normal, repulsConst, screenLength);
        push_back(b);
        return b ;
    }
    
    ///create a spherical boundary element
    BoundaryElement* createSphereBoundaryElement(vector<double>& coords, double radius, double repulsConst, double screenLength) {
        
        BoundaryElement* b = new SphereBoundaryElement(coords, radius, repulsConst, screenLength);
        push_back(b);
        return b ;
    }
    
    // Remove boundary element
    void removeBoundaryElement(BoundaryElement* b){
        
        remove(b);
        delete b;
        
    }
private:
    static BoundaryElementDB* _instance;
    BoundaryElementDB() {};
    
};








#endif /* defined(__Cyto__BoundaryElementDB__) */
