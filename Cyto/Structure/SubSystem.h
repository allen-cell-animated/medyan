//
//  SubSystem.h
//  CytoMech
//
//  Created by Konstantin Popov on 4/15/14.
//  Copyright (c) 2014 Konstantin Popov. All rights reserved.
//

#ifndef __CytoMech__SubSystem__
#define __CytoMech__SubSystem__

#include <iostream>
#include <vector>
#include <unordered_set>

#include "common.h"

///FORWARD DECLARATIONS
class Boundary;
class Filament;
class Cylinder;
class Linker;
class MotorGhost;

class Movable;
class Reactable;

/*! This is the main class which handles all changes and information regarding the system. 
    This class operates as a top manager and provides connections between smaller parts of the system. 
    All parts crations and chenges go through this class and will be redirected to lower levels.
 */
class SubSystem {
public:

    ///Add a boundary to this subsystem
    void addBoundary(Boundary* boundary) {_boundary = boundary;}
    
    /// Add new Filaments. v - coordinates of the first and last bead in the filament.
    void addNewFilaments(vector<vector<vector<double>>>& v);
    /// Add a new filament at runtime
    void addNewFilament(vector<vector<double>>& v);
    ///remove a filament from the system
    void removeFilament(Filament* f);
    
    /// Add a linker at initialization
    void addNewLinkers(vector<vector<Cylinder*>> &v, short linkerType);
    /// Add a single linker during runtime
    void addNewLinker(Cylinder* c1, Cylinder* c2, short linkerType, double position1, double position2);
    /// remove a linker from the system
    void removeLinker(Linker* l);
    
    //Add many motors at initialization
    void addNewMotorGhosts(vector<vector<Cylinder*>>& v, short motorType);
    /// Add a motor during runtime
    void addNewMotorGhost(Cylinder* c1, Cylinder* c2, short motorType, double position1, double position2);
    /// remove a motor ghost from the system
    void removeMotorGhost(MotorGhost* m);
    
    ///add and remove movables
    void addMovable(Movable* mov) { _movables.insert(mov); }
    void removeMovable(Movable* mov) {
        auto it = _movables.find(mov);
        if(it != _movables.end()) _movables.erase(it);
    }
    const unordered_set<Movable*>& getMovables() {return _movables;}
    
    ///add and remove reactables
    void addReactable(Reactable* r) { _reactables.insert(r); }
    void removeReactable(Reactable* r) {
        auto it = _reactables.find(r);
        if(it != _reactables.end()) _reactables.erase(it);
    }
    const unordered_set<Reactable*>& getReactables() {return _reactables;}
    
    //System related iterfaces:
    int getSystemSize(); //Return a number of beads;
    double getSubSystemEnergy(); // Return a value of the parameter _energy (NOT COMPUTIN!).
    
    Boundary* getBoundary() {return _boundary;}
    
    void setSubSystemEnergy(double energy); //set energy of subsystem
private:
    double _energy = 0; ///< energy of subsystem
    Boundary* _boundary; ///<subsystem boundary
    
    unordered_set<Movable*> _movables;
    unordered_set<Reactable*> _reactables;
};



#endif /* defined(__CytoMech__SubSystem__) */
