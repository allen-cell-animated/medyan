//
//  NeighborList.h
//  Cyto
//
//  Created by James Komianos on 11/12/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#ifndef __Cyto__NeighborList__
#define __Cyto__NeighborList__

#include <stdio.h>

#include <stdio.h>
#include <unordered_map>
#include <vector>

#include "common.h"

///FORWARD DECLARATIONS
class Neighbor;

class NeighborList {
    
protected:
    unordered_map<Neighbor*, vector<Neighbor*>> _list; ///The neighbors list, as a hash map
    float _rMax;  ///< max distance cutoff
    float _rMin;  ///< min distance cutoff
    
public:
    ///Constructor and destructor
    NeighborList(float rMax, float rMin) : _rMax(rMax), _rMin(rMin) {}
    
    ///Destructor
    /// @note noexcept is important here. Otherwise, gcc flags the constructor as potentially throwing,
    /// which in turn disables move operations by the STL containers. This behaviour is a gcc bug
    /// (as of gcc 4.703), and will presumbaly be fixed in the future.
    virtual ~NeighborList() noexcept {}
    
    ///re-initialize the neighborlist
    virtual void reset() {
        //loop through all neighbor keys
        for(auto it = _list.begin(); it != _list.end(); it++) {
            
            it->second.clear(); ///clear vector of neighbors
            updateNeighbors(it->first);
        }
    }
    
    ///add and remove a neighbor
    virtual void removeNeighbor(Neighbor* n) {_list.erase(n);}
    
    ///get all neighbors of a given
    virtual const vector<Neighbor*>& getNeighbors(Neighbor* n) {return _list[n];}
    
    ///Add and update neighbors
    virtual void addNeighbor(Neighbor* n) = 0;
    virtual void updateNeighbors(Neighbor* n) = 0;
    
};

class CylinderNeighborList : public NeighborList {
    
private:
    bool _crossFilamentOnly; ///< whether to include cylinders in same filament
    
public:
    ///Constructor and destructor
    CylinderNeighborList(float rMax, float rMin, bool crossFilamentOnly = false) :
                         NeighborList(rMax, rMin), _crossFilamentOnly(crossFilamentOnly) {}
    ~CylinderNeighborList() {}
    
    ///add and remove a neighbor
    virtual void addNeighbor(Neighbor* n);
    
    ///Update neighbors of a given neighbor
    virtual void updateNeighbors(Neighbor* n);

    
};

class BoundaryElementNeighborList : public NeighborList {
    
public:
    ///Constructor and destructor
    BoundaryElementNeighborList(float rMax, float rMin) : NeighborList(rMax, rMin){}
    ~BoundaryElementNeighborList() {}
    
    ///add and remove a neighbor
    virtual void addNeighbor(Neighbor* n);
    
    ///Update neighbors of a given neighbor
    virtual void updateNeighbors(Neighbor* n);
};



#endif /* defined(__Cyto__NeighborList__) */
