//
//  Filament.h
//  CytoMech
//
//  Created by Konstantin Popov on 4/15/14.
//  Copyright (c) 2014 Konstantin Popov. All rights reserved.
//

#ifndef __CytoMech__Filament__
#define __CytoMech__Filament__

#include <vector>
#include <iostream>
#include <deque>
#include <math.h>
#include "Cylinder.h"

class SubSystem;
///Filament class is used to store data about connectivity of cylinders
/*!
 * This class contains information about beads connectivity. It iterates over the beads whithin itself
 * and compute a bonded contribution to forces and energy of beads.
 */


class Filament {

private:
    std::deque<Cylinder*> _pCylinderVector; //< Vector of cylinders;
    Cylinder* _pLastCylinder = nullptr;
    SubSystem* _pSubSystem;
    
public:
    /// This constructor creates a short filament, containing only two beads. Coordinates of the first bead is an
    /// input, second is set up by using an input direction and a coarsegrain length L. Using all this, two constructors
    /// for beads are called.
	Filament(SubSystem* ps, std::vector<double> position, std::vector<double> direction);
    /// This constructor is called to create a longer filament. It creates a filament with a number of beads numBeads.
    /// Filaments starts and ends in the point determined by position vector and has a direction direction. Number of
    /// beads is equal to the number of cylinders. The last cylinder doesnt have an end(second) bead and will not be
    /// pushed to cylinder vector, but will be stored in the _pLastCylinder;
    Filament(SubSystem* ps, std::vector<std::vector<double>> position, int numBeads);
    
    void PolymerizeFront();  // Polymerization of a new cylinder (with the first bead = new bead b and last bead = empty: before: --x---x---o, after --x---x---[x---o]). Next position is based on previous beads directions in the filament. This function creates a new bead. So, this function mostly called during further polymerization, not initiation.
    void PolymerizeBack();  // Same as Polymerization frond, but adds a new first cylinder with first bead = new bead and a second bead is equal to the firs bead in the cylynder, which used to be first.
    
    ///Polymerize, used for initialization
    void PolymerizeFront(std::vector<double> coordonates );
    void PolymerizeBack(std::vector<double> coordonates );
    
    ///Delete a bead from this filament
    void DeleteBead(Bead*);
    
    ///More getters
    Cylinder* getLastCylinder() {return _pLastCylinder;}
    
    std::deque<Cylinder*> getCylinderVector() {return _pCylinderVector;}
    
    //just for example, will rewrite this function, so it not returns anything
    std:: vector<double> NextBeadProjection(Bead* pb, double d, std::vector<double> director);
    std::vector<std::vector<double> > StraightFilamentProjection(std::vector<std::vector<double>> v, int numBeads);
    std:: vector<std::vector<double> > ArcFilamentProjection(std::vector<std::vector<double>> v, int numBeads);
    
};

#endif /* defined(__CytoMech__Filament__) */
