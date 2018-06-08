
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v3.1
//
//  Copyright (2015-2016)  Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the MEDYAN web page for more information:
//  http://www.medyan.org
//------------------------------------------------------------------

#ifndef MEDYAN_NeighborListImpl_h
#define MEDYAN_NeighborListImpl_h

#include <unordered_map>
#include <vector>

#include "common.h"

#include "NeighborList.h"
#include "DynamicNeighbor.h"
#include "BinGrid.h"
#include "SysParams.h"

//FORWARD DECLARATIONS
class Cylinder;
class Bubble;
class BoundaryElement;

/// An implementation of NeighborList for Cylinder-Cylinder interactions
/// This can be a half or full list depending on the usage.
class CylinderCylinderNL : public NeighborList {
    
private:
    unordered_map<Cylinder*, vector<Cylinder*>> _list;
    ///< The neighbors list, as a hash map
    
    bool _full; ///<Specifying whether this is a full or half list

#ifdef CUDAACCL_NL
    vector<int> blocksnthreads;
    int nint;
//    vector<int> pair_cIndex_cmp;
    vector<int> pair_cIndex_cnIndex;
    int *gpu_pair_cIndex_cnIndex;
//    int *gpu_pair_cIndex_cmp;
    double *gpu_params = NULL;
    int *gpu_NL = NULL;
    int *gpu_numpairs = NULL;
    int *gpu_params2;
    int numpairs[1];
#endif
    
    ///Helper function to update neighbors
    ///@param runtime - specifying whether the cylinder is being
    ///created/destroyed at runtime vs at a full neighbor list update.
    void updateNeighbors(Cylinder* cylinder, bool runtime = false);

    
public:
#ifdef CUDAACCL_NL
    bool cudacpyforces = false;

    int getNLsize() {
        return numpairs[0];
    }
    int* getNLsizeCUDA(){
        return gpu_numpairs;
    }
    int* getNLCUDA(){
        return gpu_NL;
    }
#endif
    short _ID; //ID helps link binGridType to NeighborList.

    CylinderCylinderNL(float rMax, float rMin = 0.0, bool full = false, short ID = 0)
            : NeighborList(rMax, rMin), _full(full) {
    }
    virtual void addNeighbor(Neighbor* n);
    virtual void removeNeighbor(Neighbor* n);
    
    //@{
    /// The implementation of these functions calls the static version,
    /// all cylinders are dynamic
    virtual void addDynamicNeighbor(DynamicNeighbor* n) {addNeighbor(n);}
    virtual void removeDynamicNeighbor(DynamicNeighbor* n) {removeNeighbor(n);}
    //@}
    
    virtual void reset();
    
    /// Get all cylinder neighbors
    vector<Cylinder*> getNeighbors(Cylinder* cylinder);
    
};

/// An implementation of NeighborList for BoundaryElement-Cylinder interactions
class BoundaryCylinderNL : public NeighborList {
    
private:
    unordered_map<BoundaryElement*, vector<Cylinder*>> _list;
    ///< The neighbors list, as a hash map
    
    ///Helper function to update neighbors
    void updateNeighbors(BoundaryElement* be);
    
public:
    BoundaryCylinderNL(float rMax): NeighborList(rMax) {}
    
    virtual void addNeighbor(Neighbor* n);
    virtual void removeNeighbor(Neighbor* n);
    
    virtual void addDynamicNeighbor(DynamicNeighbor* n);
    virtual void removeDynamicNeighbor(DynamicNeighbor* n);
    
    virtual void reset();
    
    /// Get all Cylinder neighbors of a boundary element
    vector<Cylinder*> getNeighbors(BoundaryElement* be);
};


/// An implementation of NeighborList for BoundaryElement-Bubble interactions
class BoundaryBubbleNL : public NeighborList {
    
private:
    unordered_map<BoundaryElement*, vector<Bubble*>> _list;
    ///< The neighbors list, as a hash map
    
    ///Helper function to update neighbors
    void updateNeighbors(BoundaryElement* be);
    
public:
    BoundaryBubbleNL(float rMax): NeighborList(rMax) {}
    
    virtual void addNeighbor(Neighbor* n);
    virtual void removeNeighbor(Neighbor* n);
    
    virtual void addDynamicNeighbor(DynamicNeighbor* n);
    virtual void removeDynamicNeighbor(DynamicNeighbor* n);
    
    virtual void reset();
    
    /// Get all Bubble neighbors of a boundary element
    vector<Bubble*> getNeighbors(BoundaryElement* be);
};

/// An implementation of NeighborList for Bubble-Bubble interactions
/// @note - This is currently implemented as a half list only
class BubbleBubbleNL : public NeighborList {
    
private:
    unordered_map<Bubble*, vector<Bubble*>> _list;
    ///< The neighbors list, as a hash map
    
    ///Helper function to update neighbors
    void updateNeighbors(Bubble* bb);
    
public:
    BubbleBubbleNL(float rMax): NeighborList(rMax) {}
    
    virtual void addNeighbor(Neighbor* n);
    virtual void removeNeighbor(Neighbor* n);
    
    //@{
    /// The implementation of these functions calls the static version,
    /// all Bubbles are dynamic
    virtual void addDynamicNeighbor(DynamicNeighbor* n) {addNeighbor(n);}
    virtual void removeDynamicNeighbor(DynamicNeighbor* n) {removeNeighbor(n);}
    //@}
    
    virtual void reset();
    
    /// Get all Bubble neighbors of a bubble
    vector<Bubble*> getNeighbors(Bubble* bb);
};

/// An implementation of NeighborList for Bubble-Cylinder interactions
class BubbleCylinderNL : public NeighborList {
    
private:
    unordered_map<Bubble*, vector<Cylinder*>> _list;
    ///< The neighbors list, as a hash map
    
    ///Helper function to update neighbors
    void updateNeighbors(Bubble* bb);
    
public:
    BubbleCylinderNL(float rMax): NeighborList(rMax) {}
    
    virtual void addNeighbor(Neighbor* n);
    virtual void removeNeighbor(Neighbor* n);
    
    //@{
    /// The implementation of these functions calls the static version,
    /// all Bubbles and Cylinders are dynamic
    virtual void addDynamicNeighbor(DynamicNeighbor* n) {addNeighbor(n);}
    virtual void removeDynamicNeighbor(DynamicNeighbor* n) {removeNeighbor(n);}
    //@}
    
    virtual void reset();
    
    /// Get all Cylinder neighbors of a bubble
    vector<Cylinder*> getNeighbors(Bubble* bb);
};





#endif
