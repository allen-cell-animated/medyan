//
//  NeighborsManager.h
//  Cyto
//
//  Created by James Komianos on 11/10/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#ifndef __Cyto__NeighborsManager__
#define __Cyto__NeighborsManager__

#include <stdio.h>

#include "common.h"
#include "MathFunctions.h"

#include "CylinderDB.h"
#include "GController.h"

using namespace mathfunc;

///FORWARD DECLARATIONS
class Cylinder;


class NeighborList {
    
private:
    unordered_map<Cylinder*, vector<Cylinder*>> _list; ///The neighbors list, as a hash map
    float _rMax;
    float _rMin;
    
    vector<Cylinder*> findNearbyCylinders(Cylinder* cylinder);
    
public:
    
    NeighborList(float rMax, float rMin) :  _rMax(rMax), _rMin(rMin) { reset(); }
    ~NeighborList() {}
    
    void reset();
    
};




class NeighborListManager {
    
private:
    vector<NeighborList> _chemNeighborLists;
    vector<NeighborList> _mechNeighborLists;

public:
    
    NeighborListManager() {}
    ~NeighborListManager() {}
    

    void addChemNeighborList(float rMax = 0.0, float rMin = 0.0) {
        
        _chemNeighborLists.emplace_back(rMax, rMin);
    }
    
    void addMechNeighborList(float rMax = 0.0, float rMin = 0.0) {
        
        _mechNeighborLists.emplace_back(rMax, rMin);
    }
    
    
    

    
    
    
    
    
    
    
    
};



#endif /* defined(__Cyto__NeighborsManager__) */
