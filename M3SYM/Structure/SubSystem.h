 
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

#include "NeighborList.h"
#include "DynamicNeighbor.h"

#include "SysParams.h"

//FORWARD DECLARATIONS
class Boundary;
class Filament;
class Cylinder;
class Linker;
class MotorGhost;
class BranchingPoint;

class CompartmentGrid;

/// Manages all [Movables](@ref Movable) and [Reactables](@ref Reactable). Also holds all
/// [NeighborLists](@ref NeighborList) associated with chemical or mechanical interactions,
/// as well as the CompartmentGrid which contains all chemical structural information, and the
/// system Boundary.

/*! This is a class which handles all changes and information regarding the simulated system.
 *  This class operates as a top manager and provides connections between smaller parts 
 *  of the system. All creation and changes go through this class and will be redirected 
 *  to lower levels. See the databases for more documentation on the explicit creation of
 *  subsystem objects at initialization and during runtime.
 *
 *  This class has functions to add or remove Trackable elements from the system, as well
 *  as update Movable and Reactable instances in the system. It also can update the 
 *  NeighborList container that it holds for the system.
 */
class SubSystem {
    
public:
    ///Default constructor does nothing
    SubSystem() {} ~SubSystem() {}
    
    /// Add a Trackable to the SubSystem
    template<class T, typename ...Args>
    T* addTrackable(Args&& ...args) {
        
        //create instance
        T* t = new T( forward<Args>(args)...); t->addToSubSystem();
        
        //if movable or reactable, add
        if(t->_movable) addMovable((Movable*)t);
        
        if(t->_reactable) addReactable((Reactable*)t);
        
        //if neighbor, add
        if(t->_dneighbor) {
            for(auto nlist : _neighborLists.getElements())
                nlist->addDynamicNeighbor((DynamicNeighbor*)t);
        }
        
        else if(t->_neighbor) {
            for(auto nlist : _neighborLists.getElements())
                nlist->addNeighbor((Neighbor*)t);
        }
        
        return t;
    }
    
    /// Remove a trackable from the SubSystem
    template<class T>
    void removeTrackable(T* t) {
        
        //remove from subsystem
        t->removeFromSubSystem();
        
        //if movable or reactable, remove
        if(t->_movable) removeMovable((Movable*)t);
        
        if(t->_reactable) removeReactable((Reactable*)t);
        
        //if neighbor, remove
        if(t->_dneighbor) {
            for(auto nlist : _neighborLists.getElements())
                nlist->removeDynamicNeighbor((DynamicNeighbor*)t);
        }
        
        else if(t->_neighbor) {
            for(auto nlist : _neighborLists.getElements())
                nlist->removeNeighbor((Neighbor*)t);
        }
        
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
  
    /// Reset all neighbor lists in subsystem
    void resetNeighborLists() {
        
        for(auto nl: _neighborLists.getElements())
            nl->reset();
    }
    
    //@{
    ///Subsystem energy management
    double getSubSystemEnergy() {return _energy;}
    void setSubSystemEnergy(double energy) {_energy = energy;}
    //@}
    
    //@{
    /// CompartmentGrid management
    void setCompartmentGrid(CompartmentGrid* grid) {_compartmentGrid = grid;}
    CompartmentGrid* getCompartmentGrid() {return _compartmentGrid;}
    //@]
    
    /// Update the binding managers of the system
    void updateBindingManagers();
    
private:
    double _energy = 0; ///< Energy of this subsystem
    Boundary* _boundary; ///< Boundary pointer
    
    unordered_set<Movable*> _movables; ///< All movables in the subsystem
    unordered_set<Reactable*> _reactables; ///< All reactables in the subsystem
        
    Database<NeighborList*> _neighborLists; ///< All neighborlists in the system
        
    CompartmentGrid* _compartmentGrid; ///< The compartment grid
};

#endif
