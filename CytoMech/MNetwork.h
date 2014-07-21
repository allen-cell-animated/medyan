//
//  MNetwork.h
//  CytoMech
//
//  Created by Konstantin Popov on 4/15/14.
//  Copyright (c) 2014 Konstantin Popov. All rights reserved.
//

#ifndef __CytoMech__MNetwork__
#define __CytoMech__MNetwork__

#include <iostream>
#include <list>
#include <vector>
#include "MFilamentDB.h"
#include "MLinkerDB.h"
#include "MMotorGhostDB.h"



class Network
{
    
public:
    
    Network(System* ps); //!< Constructor.
    
    void ResetForces(BeadDB&); //!< Set force to zero before every update.
    void ResetForcesAux(BeadDB&); //!< Set forceAux to zero before every update.
    void CopmuteForce(int);         //!< compute forces (also an interface whic calls compute at the network sub units). Argument determins where to write forces: foce (0), forceAux (other then 0).
    double CopmuteEnergy(double);   //!< compute energy (also an interface whic calls compute at the network sub units). Argument determins if use just energy equations with a simple coordinates argument(d=0.0) or Bead.coordinates - d*Bead.force as an argument.

// Interfaces to add new objects:
    
    
    void AddNewFilaments(std::vector<std::vector<std::vector<double> > > v);    //!< Add new Filaments. v - coordinates of the first and last bead in the filament.
    void AddNewLinkers(std::vector<std::vector<Cylinder* > > v, double stretchConst);   //!< Add a linker conecting two beads Bead*[0] and Bead*[1] and providing a stretching spring constant(coubling to twisting will be implemented later).
    void AddNewLinker(Cylinder* pc1, Cylinder* pc2, double stretchConst);
    
    void AddNewMotorGhost(Cylinder* pc1, Cylinder* pc2, double k, double position1, double position2); //!< Add a ghost motor connectind two CG segments (4 beads), providing a stretching constant, and positions on each segment.
    void AddNewMotorGhosts(std::vector<std::vector<Cylinder* > > v, double k, double position1, double position2); //Add many motors. Input: vector of vectors(cyl1, cyl2), pair connected by a motor. 
    
    
    FilamentDB* getFilamentDB() {return _pfdb;}
private:
    
    System* _pS;
    FilamentDB* _pfdb;
    LinkerDB* _pldb;
    MotorGhostDB* _pmgdb;
    
    
};


#endif /* defined(__CytoMech__MNetwork__) */
