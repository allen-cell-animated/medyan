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

class Boundary;
class Cylinder;
class Compartment;

/*! This is the main class which handles all changes and information regarding the system. This class operates as a top manager and provides connections between smaller parts of the system. All parts crations and chenges go through this class and will be redirected to lower levels.
 */
class SubSystem {
public:
    // Interfaces to add new objects:
    
    ///Add a boundary to this subsystem
    void AddBoundary(Boundary* boundary) {_boundary = boundary;}
    
    /// Add new Filaments. v - coordinates of the first and last bead in the filament.
    void AddNewFilaments(std::vector<std::vector<std::vector<double>>>& v);
    
    /// Add a linker at initialization
    void AddNewLinkers(std::vector<std::vector<Cylinder*>> &v, short linkerType);
    /// Add a single linker during runtime
    void AddNewLinker(Cylinder* pc1, Cylinder* pc2, short linkerType, double position1, double position2);
    
    //Add many motors at initialization
    void AddNewMotorGhosts(std::vector<std::vector<Cylinder*>>& v, short motorType);
    /// Add a motor during runtime
    void AddNewMotorGhost(Cylinder* pc1, Cylinder* pc2, short motorType, double position1, double position2);
    
    //System related iterfaces:
    int getSystemSize(); //Return a number of beads;
    double getSubSystemEnergy(); // Return a value of the parameter _energy (NOT COMPUTIN!).
    
    void setSubSystemEnergy(double energy); //set energy of subsystem
private:
    double _energy = 0; ///< energy of subsystem
    Boundary* _boundary; ///<subsystem boundary
	
};



#endif /* defined(__CytoMech__SubSystem__) */
