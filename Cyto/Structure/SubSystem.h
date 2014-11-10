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
#include <list>
#include <vector>

#include "common.h"

#include "FilamentDB.h"
#include "LinkerDB.h"
#include "MotorGhostDB.h"

class Boundary;
class Cylinder;
class Compartment;

/*! This is the main class which handles all changes and information regarding the system. This class operates as a top manager and provides connections between smaller parts of the system. All parts crations and chenges go through this class and will be redirected to lower levels.
 */
class SubSystem {
public:
    // Interfaces to add new objects:
    
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
    void addNewLinker(Cylinder* pc1, Cylinder* pc2, short linkerType, double position1, double position2);
    /// remove a linker from the system
    void removeLinker(Linker* l);
    
    //Add many motors at initialization
    void addNewMotorGhosts(vector<vector<Cylinder*>>& v, short motorType);
    /// Add a motor during runtime
    void addNewMotorGhost(Cylinder* pc1, Cylinder* pc2, short motorType, double position1, double position2);
    /// remove a motor ghost from the system
    void removeMotorGhost(MotorGhost* m);
    
    
    //System related iterfaces:
    int getSystemSize(); //Return a number of beads;
    double getSubSystemEnergy(); // Return a value of the parameter _energy (NOT COMPUTIN!).
    
    Boundary* getBoundary() {return _boundary;}
    
    void setSubSystemEnergy(double energy); //set energy of subsystem
private:
    double _energy = 0; ///< energy of subsystem
    Boundary* _boundary; ///<subsystem boundary
	
};



#endif /* defined(__CytoMech__SubSystem__) */
