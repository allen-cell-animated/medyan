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
#include "Neighbor.h"

///FORWARD DECLARATIONS
class Cylinder;

class NeighborList {
    
protected:
    unordered_map<Neighbor*, vector<Neighbor*>> _list; ///The neighbors list, as a hash map
    
    float _rMax;  ///< max distance cutoff
    float _rMin;  ///< min distance cutoff
    
    short _ID; ///< unique ID of this neighborlist (corresponding to a rxn or force field)
    
public:
    ///Constructor and destructor
    NeighborList(float rMax, float rMin, int ID) : _rMax(rMax), _rMin(rMin), _ID(ID) {}
    
    ///Destructor
    /// @note noexcept is important here. Otherwise, gcc flags the constructor as potentially throwing,
    /// which in turn disables move operations by the STL containers. This behaviour is a gcc bug
    /// (as of gcc 4.703), and will presumbaly be fixed in the future.
    virtual ~NeighborList() noexcept {}
    
    ///re-initialize the neighborlist
    virtual void reset() = 0;
    
    virtual void addNeighbor(Neighbor* n) = 0;
    virtual void removeNeighbor(Neighbor* n) = 0;
    
    const unordered_map<Neighbor*, vector<Neighbor*>>& getList() {return _list;}
    
};

class CylinderNeighborList : public NeighborList {
    
private:
    bool _crossFilamentOnly; ///< whether to include cylinders in same filament
    
    ///Find all nearby cylinders relative to a given cylinder
    vector<Cylinder*> findNearbyCylinders(Cylinder* cylinder);
    
public:
    ///Constructor and destructor
    CylinderNeighborList(float rMax, float rMin, int ID, bool crossFilamentOnly = false) :
                         NeighborList(rMax, rMin, ID), _crossFilamentOnly(crossFilamentOnly) { reset(); }
    ~CylinderNeighborList() {}
    
    ///Re-initialize a neighbors list
    virtual void reset();
    
    virtual void addNeighbor(Neighbor* n);
    virtual void removeNeighbor(Neighbor* n);
    
};


#endif /* defined(__Cyto__NeighborList__) */
