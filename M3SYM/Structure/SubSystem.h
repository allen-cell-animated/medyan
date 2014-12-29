
//------------------------------------------------------------------
//  **M3SYM** - Simulation Package for the Mechanochemical
//              Dynamics of Active Networks, 3rd Generation
//
//  Copyright (2014) Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the Papoian lab page for installation and documentation:
//  http://papoian.chem.umd.edu/
//------------------------------------------------------------------

#ifndef M3SYM_SubSystem_h
#define M3SYM_SubSystem_h

#include <vector>
#include <unordered_set>

#include "common.h"
#include "NeighborListContainer.h"

#include "SystemParameters.h"

//FORWARD DECLARATIONS
class Boundary;
class Filament;
class Cylinder;
class Linker;
class MotorGhost;
class BranchingPoint;

class Movable;
class Reactable;


/// Manages all objects in the system, including [Filaments] (@ref Filament), [Linkers]
/// (@ref Linker), [MotorGhosts] (@ref MotorGhost), and [BranchingPoints](@ref
/// BranchingPoint).

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
    SubSystem() : CBENLContainer(SystemParameters::Boundaries().boundaryCutoff) {}
#endif
        
    /// Add new [Filaments](@ref Filament).
    /// @param v - coordinates of the first and last bead in the filament.
    void addNewFilaments(vector<vector<vector<double>>>& v);
    /// Add a new Filament at runtime
    Filament* addNewFilament(vector<vector<double>>& v);
    /// Remove a Filament from the system
    void removeFilament(Filament* f);
    
    /// Add [Linkers](@ref Linker) at initialization
    /// @param v - vector of cylinders to connect to
    void addNewLinkers(vector<vector<Cylinder*>> &v, short linkerType);
    /// Add a single Linker during runtime
    Linker* addNewLinker(Cylinder* c1, Cylinder* c2, short linkerType,
                         double position1, double position2);
    /// Remove a Linker from the system
    void removeLinker(Linker* l);
    
    /// Add [MotorGhosts](@ref MotorGhost) at initialization
    /// @param v - vector of cylinders to connect to
    void addNewMotorGhosts(vector<vector<Cylinder*>>& v, short motorType);
    /// Add a MotorGhost during runtime
    MotorGhost* addNewMotorGhost(Cylinder* c1, Cylinder* c2, short motorType,
                                 double position1, double position2);
    /// remove a MotorGhost ghost from the system
    void removeMotorGhost(MotorGhost* m);
    
    /// Add [BranchingPoints](@ref BranchingPoint) at initialization
    /// @param v - vector of cylinders to connect to
    void addNewBranchingPoints(vector<vector<Cylinder*>>& v, short branchType);
    /// Add a BranchingPoint during runtime
    BranchingPoint* addNewBranchingPoint(Cylinder* c1, Cylinder* c2,
                                         short branchType, double position);
    /// remove a BranchingPoint from the system
    void removeBranchingPoint(BranchingPoint* b);
    
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
    
    /// Return the number of [Beads](@ref Bead) in the system
    int getSystemSize();
    
    //@{
    /// Subsystem energy management
    double getSubSystemEnergy();
    void setSubSystemEnergy(double energy);
    //@}
    
    /// Get the subsystem boundary
    Boundary* getBoundary() {return _boundary;}
    /// Add a boundary to this subsystem
    void addBoundary(Boundary* boundary) {_boundary = boundary;}
    
    /// Get the cylinders that are currently interacting with a boundary
    vector<Cylinder*> getBoundaryCylinders();
    
private:
    double _energy = 0; ///< energy
    Boundary* _boundary; ///< boundary pointer
    
    unordered_set<Movable*> _movables; ///< All movables in the subsystem
    unordered_set<Reactable*> _reactables; ///< All reactables in the subsystem
};

#endif
