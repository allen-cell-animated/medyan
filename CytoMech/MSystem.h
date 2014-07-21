//
//  MSystem.h
//  CytoMech
//
//  Created by Konstantin Popov on 4/15/14.
//  Copyright (c) 2014 Konstantin Popov. All rights reserved.
//

#ifndef __CytoMech__MSystem__
#define __CytoMech__MSystem__

#include <iostream>
#include <list>
#include <vector>
#include "MBeadDB.h"
#include "MBead.h"
#include "MCylinder.h"
#include "MCylinderDB.h"
#include "MFilamentDB.h"
#include "Mcommon.h"

class System
{
    
    /*! This is the main class which handel all changes and information regarding the system. This class operates as a top manager and provides connections between smaller parts of the system. All parts crations and chenges go through this class and will be redirected to lower levels.
     
     */
    
public:
	
    System();
	
    //"Creational" interfaces:
	Bead* AddNewBead(std::vector<double> v); //Add a new bead with coordinates v ;
	
	Bead* AddNewBead(); // Add an empty bead. A dummy interface;
	
    Cylinder* AddNewCylinder(Filament* pf, Bead* pb); // Create a cylinder on the filament pf with the first bead pb;
    
    //	void AddNewFilament(vector<double> v); //A dummmy interface;
    
    void AddNewFilaments(std::vector<std::vector<std::vector<double> > > v, int i); // An interface that adds new filaments, where vector v is a multidim vector with coordinates of the first and last bead of the filament, and int i is a number of network where to create it (in case of multiple networks).
    
    Network* AddNewNetwork(); // Add a new network. When Sysytem is created it automatically creates one network;
    
    void AddNewLinkers(std::vector<std::vector<Cylinder* > > v, double stretchConst, int i);
    void AddNewLinker(Cylinder* pc1, Cylinder* pc2, double stretchConst, int i);
    
    void AddNewMotorGost(Cylinder* pc1, Cylinder* pc2, double k, double position1, double position2, int i);
    void AddNewMotorGhosts(std::vector<std::vector<Cylinder* > > v, double k, double position1, double position2, int i);
    
    //Destructional iterfaces:
    
    void RemoveBead(Bead*); // Not fully impemented yet.
    void SplitFilament(Filament*, Bead*); // Not fully implimented yet.
    
    //System related iterfaces:
    
	BeadDB* getBDB();   //Return a pointer to a BeadDB;
    CylinderDB* getCDB(); //Return a pointer to a CylinderDB;
    int getSystemSize(); //Return a number of beads;
    double getSystemEnergy(); // Return a value of the parameter _energy (NOT COMPUTIN!).
    
    std::vector<Network*> getNetwork() {return _pNetworkVector;}
    
    
    // Mechanical interfaces:
    void CopmuteForce(int); //Interface which goes over all sub parts of the system and starts force calculations in them. They results will be written in Bead.force based on Bead.coordinates.
    double UpdateEnergy(double); //Compute and update a total energy.
    
    
    
private:
	BeadDB* _pbdb;
    CylinderDB* _pcdb;
    double _energy;
	
    std::vector<Network*> _pNetworkVector;
	// return filament vector here!!!
};



#endif /* defined(__CytoMech__MSystem__) */
