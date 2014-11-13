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
class Cylinder;
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
    virtual void reset() = 0;
    
    ///add and remove a neighbor
    virtual void addNeighbor(Neighbor* n) = 0;
    virtual void removeNeighbor(Neighbor* n) = 0;
    
    ///Update neighbors
    virtual void updateNeighbors(Neighbor* n) = 0;
    
    ///get all neighbors of a given
    virtual const vector<Neighbor*>& getNeighbors(Neighbor* n) = 0;
    
};

class CylinderNeighborList : public NeighborList {
    
private:
    bool _crossFilamentOnly; ///< whether to include cylinders in same filament
    
public:
    ///Constructor and destructor
    CylinderNeighborList(float rMax, float rMin, bool crossFilamentOnly = false) :
                         NeighborList(rMax, rMin), _crossFilamentOnly(crossFilamentOnly) {}
    ~CylinderNeighborList() {}
    
    ///Re-initialize a neighbors list
    virtual void reset();
    
    ///add and remove a neighbor
    virtual void addNeighbor(Neighbor* n);
    virtual void removeNeighbor(Neighbor* n);
    
    ///Update neighbors of a given neighbor
    virtual void updateNeighbors(Neighbor* n);
    
    virtual const vector<Neighbor*>& getNeighbors(Neighbor* n);
    
};


#endif /* defined(__Cyto__NeighborList__) */
