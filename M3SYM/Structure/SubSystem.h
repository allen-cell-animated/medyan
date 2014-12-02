
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

//FORWARD DECLARATIONS
class Boundary;
class Filament;
class Cylinder;
class Linker;
class MotorGhost;

class Movable;
class Reactable;


/// Manages all objects in the system, including [Filaments] (@ref Filament), [Linkers] (@ref Linker), and [Motors] (@ref Motor).

/*! This is a class which handles all changes and information regarding the system.
 *  This class operates as a top manager and provides connections between smaller parts of the system.
 *  All creation and changes go through this class and will be redirected to lower levels. See databases for more
 *  documentation on the explicit creation of subsystem objects at initialization and during runtime.
 */
class SubSystem {
public:

    /// Add a boundary to this subsystem
    void addBoundary(Boundary* boundary) {_boundary = boundary;}
    
    /// Add new filaments.
    /// @param v - coordinates of the first and last bead in the filament.
    void addNewFilaments(vector<vector<vector<double>>>& v);
    /// Add a new filament at runtime
    void addNewFilament(vector<vector<double>>& v);
    /// Remove a filament from the system
    void removeFilament(Filament* f);
    
    /// Add a linker at initialization
    /// @param v - vector of cylinders to connect to
    void addNewLinkers(vector<vector<Cylinder*>> &v, short linkerType);
    /// Add a single linker during runtime
    void addNewLinker(Cylinder* c1, Cylinder* c2, short linkerType, double position1, double position2);
    /// Remove a linker from the system
    void removeLinker(Linker* l);
    
    /// Add a motor at initialization
    /// @param v - vector of cylinders to connect to
    void addNewMotorGhosts(vector<vector<Cylinder*>>& v, short motorType);
    /// Add a motor during runtime
    void addNewMotorGhost(Cylinder* c1, Cylinder* c2, short motorType, double position1, double position2);
    /// remove a motor ghost from the system
    void removeMotorGhost(MotorGhost* m);
    
    //@{
    /// Setter functions for movables
    void addMovable(Movable* mov) { _movables.insert(mov); }
    void removeMovable(Movable* mov) {
        auto it = _movables.find(mov);
        if(it != _movables.end()) _movables.erase(it);
    }
    //@}
    /// Get all movables in the subsystem
    const unordered_set<Movable*>& getMovables() {return _movables;}
    
    //@{
    /// Setter function for reactables
    void addReactable(Reactable* r) { _reactables.insert(r); }
    void removeReactable(Reactable* r) {
        auto it = _reactables.find(r);
        if(it != _reactables.end()) _reactables.erase(it);
    }
    //@}
    /// Get all reactables in the subsystem
    const unordered_set<Reactable*>& getReactables() {return _reactables;}
    
    /// Return the number of beads in the system
    int getSystemSize();
    
    double getSubSystemEnergy();
    Boundary* getBoundary() {return _boundary;}
    
    void setSubSystemEnergy(double energy);
private:
    double _energy = 0; ///< energy of subsystem
    Boundary* _boundary; ///< subsystem boundary
    
    unordered_set<Movable*> _movables; ///< All movables in the subsystem
    unordered_set<Reactable*> _reactables; ///< All reactables in the subsystem
};



#endif
