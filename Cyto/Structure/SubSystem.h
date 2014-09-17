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
#include "BeadDB.h"
#include "CylinderDB.h"
#include "FilamentDB.h"
#include "MotorGhostDB.h"
#include "LinkerDB.h"

/*! This is the main class which handles all changes and information regarding the system. This class operates as a top manager and provides connections between smaller parts of the system. All parts crations and chenges go through this class and will be redirected to lower levels.
 */
class SubSystem {
public:
    // Interfaces to add new objects:
    
    ///Add a boundary to this subsystem
    void AddBoundary(Boundary* boundary) {_boundary = boundary;}
    
    /// Add new Filaments. v - coordinates of the first and last bead in the filament.
    void AddNewFilaments(std::vector<std::vector<std::vector<double>>>& v);
    
    /// Add a linker conecting two beads Bead*[0] and Bead*[1] and providing a stretching spring constant
    /// (coupling to twisting will be implemented later).
    void AddNewLinkers(std::vector<std::vector<Cylinder*>> &v, double stretchConst);
    /// Add a single linker
    void AddNewLinker(Cylinder* pc1, Cylinder* pc2, double stretchConst);
    
    /// Add a ghost motor connectind two CG segments (4 beads), providing a stretching constant, and positions
    /// on each segment.
    void AddNewMotorGhost(Cylinder* pc1, Cylinder* pc2, double k, double position1, double position2);
    //Add many motors. Input: vector of vectors(cyl1, cyl2), pair connected by a motor.
    void AddNewMotorGhosts(std::vector<std::vector<Cylinder*>>& v, double k, double position1, double position2);
    
    //System related iterfaces:
    int getSystemSize(); //Return a number of beads;
    double getSubSystemEnergy(); // Return a value of the parameter _energy (NOT COMPUTIN!).
    
    void setSubSystemEnergy(double energy); //set energy of subsystem
private:
    double _energy = 0; ///< energy of subsystem
    Boundary* _boundary; ///<subsystem boundary
	
};



#endif /* defined(__CytoMech__SubSystem__) */
