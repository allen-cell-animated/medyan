
//------------------------------------------------------------------
//  **M3SYM** - Simulation Package for the Mechanochemical
//              Dynamics of Active Networks, 3rd Generation
//
//  Copyright (2014) Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the Papoian lab page for more information:
//  http://papoian.chem.umd.edu/
//------------------------------------------------------------------

#ifndef M3SYM_SubSystem_h
#define M3SYM_SubSystem_h

#include <vector>
#include <unordered_set>

#include "common.h"

#include "Trackable.h"
#include "Database.h"
#include "Movable.h"
#include "Reactable.h"

#include "Neighbor.h"
#include "DynamicNeighbor.h"

#ifdef DYNAMICRATES
#include "NeighborListContainer.h"
#endif

#include "SysParams.h"

//FORWARD DECLARATIONS
class Boundary;
class Filament;
class Cylinder;
class Linker;
class MotorGhost;
class BranchingPoint;

class CompartmentGrid;

/// Manages all objects in the system, including [Filaments] (@ref Filament), [Linkers]
/// (@ref Linker), [MotorGhosts] (@ref MotorGhost), and [BranchingPoints](@ref
/// BranchingPoint). Also contains the CompartmentGrid.

/*! This is a class which handles all changes and information regarding the system.
 *  This class operates as a top manager and provides connections between smaller parts 
 *  of the system. All creation and changes go through this class and will be redirected 
 *  to lower levels. See databases for more documentation on the explicit creation of 
 *  subsystem objects at initialization and during runtime.
 *
 *  The SubSystem class also extends CBENLContainer, holding a neighbors list for
 *  [Cylinders](@ref Cylinder) near boundaries. This is used for reaction rate updating.
 */
class SubSystem
#ifdef DYNAMICRATES
    : public CBENLContainer {
#else 
    {
#endif
public:
#ifdef DYNAMICRATES 
    SubSystem() : CBENLContainer(SysParams::Boundaries().BoundaryCutoff) {}
#endif
        
    /// Add a Trackable to the SubSystem
    template<class T, typename ...Args>
    T* addTrackable(Args&& ...args) {
        
        //create instance
        T* t = new T( forward<Args>(args)...);
        
        //add to subsystem
        t->addToSubSystem();
        
        //if movable or reactable, add
        Movable* m;
        if((m = dynamic_cast<Movable*>(t))) addMovable(m);
        Reactable* r;
        if((r = dynamic_cast<Reactable*>(t))) addReactable(r);
        
        //if neighbor, add
        DynamicNeighbor* dn; Neighbor* n;
        
        if((dn = dynamic_cast<DynamicNeighbor*>(t)))
            for(auto nlist : _neighborLists.getElements())
                nlist->addDynamicNeighbor(dn);
        
        else if((n = dynamic_cast<Neighbor*>(t)))
            for(auto nlist : _neighborLists.getElements())
                nlist->addNeighbor(n);
        
        return t;
    }
        
    /// Remove a trackable from the SubSystem
    void removeTrackable(Trackable* t) {
        
        //remove from subsystem
        t->removeFromSubSystem();
        
        //if movable or reactable, remove
        Movable* m;
        if((m = dynamic_cast<Movable*>(t))) removeMovable(m);
        
        Reactable* r;
        if((r = dynamic_cast<Reactable*>(t))) removeReactable(r);
        
        //if neighbor, remove as well
        DynamicNeighbor* dn; Neighbor* n;
        
        if((dn = dynamic_cast<DynamicNeighbor*>(t)))
            for(auto nlist : _neighborLists.getElements())
                nlist->removeDynamicNeighbor(dn);
        
        else if((n = dynamic_cast<Neighbor*>(t)))
            for(auto nlist : _neighborLists.getElements())
                nlist->removeNeighbor(n);
        
        //delete it
        delete t;
    }
    
    //@{
    /// Setter functions for Movable
    void addMovable(Movable* mov) { _movables.insert(mov); }
    void removeMovable(Movable* mov) {
        auto it = _movables.find(mov);
        if(it != _movables.end()) _movables.erase(it);
    }
    //@}
    /// Get all Movable
    const unordered_set<Movable*>& getMovables() {return _movables;}
    
    //@{
    /// Setter function for Reactable
    void addReactable(Reactable* r) { _reactables.insert(r); }
    void removeReactable(Reactable* r) {
        auto it = _reactables.find(r);
        if(it != _reactables.end()) _reactables.erase(it);
    }
    //@}
    /// Get all Reactable
    const unordered_set<Reactable*>& getReactables() {return _reactables;}
    
    /// Get the subsystem boundary
    Boundary* getBoundary() {return _boundary;}
    /// Add a boundary to this subsystem
    void addBoundary(Boundary* boundary) {_boundary = boundary;}
        
    /// Add a neighbor list to the subsystem
    void addNeighborList(NeighborList* nl) {_neighborLists.addElement(nl);}
    /// Get the database of neighbor lists in subsystem
    Database<NeighborList*>& getNeighborLists() {return _neighborLists;}
        
    //@{
    ///Subsystem energy management
    double getSubSystemEnergy() {return _energy;}
    void setSubSystemEnergy(double energy) {_energy = energy;}
    //@}
        
#ifdef DYNAMICRATES
    /// Get the cylinders that are currently interacting with a boundary
    vector<Cylinder*> getBoundaryCylinders();
#endif
    
private:
    double _energy = 0; ///< energy
    Boundary* _boundary; ///< boundary pointer
    
    unordered_set<Movable*> _movables; ///< All movables in the subsystem
    unordered_set<Reactable*> _reactables; ///< All reactables in the subsystem
        
    Database<NeighborList*> _neighborLists; ///< All neighborlists in the system
        
    CompartmentGrid* _compartmentGrid; ///< The compartment grid
};

#endif
