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

#include "common.h"
#include "Cylinder.h"

///FORWARD DECLARATIONS
class SubSystem;

///Filament class is used to store data about connectivity of cylinders
/*!
 * This class contains information about beads and cylinders connectivity. It iterates over the beads whithin itself
 * and compute a bonded contribution to forces and energy of beads.
 */
class Filament {

private:
    deque<Cylinder*> _cylinderVector; ///< Vector of cylinders;
    SubSystem* _subSystem;
    
    int _ID; ///< unique integer id of this filament
    
    short _deltaPlusEnd = 0;   ///< change in filament's cylinders at plus end since last snapshot
    short _deltaMinusEnd = 0; ///< change in filament's cylinders at minus end since last snapshot
    
    short _beadIDPlusEnd = 0; ///< Bead ID to assign to next plus end bead
    short _beadIDMinusEnd = -1;  ///< Bead ID to assign to next minus end bead
    
public:
    /// This constructor creates a short filament, containing only two beads. Coordinates of the first bead is an
    /// input, second is set up by using an input direction and a coarsegrain length L. Using all this, two constructors
    /// for beads are called.
	Filament(SubSystem* s, vector<double>& position, vector<double>& direction, int ID);
    /// This constructor is called to create a longer filament. It creates a filament with a number of beads numBeads.
    /// Filaments starts and ends in the point determined by position vector and has a direction direction. Number of
    /// beads is equal to the number of cylinders. The last cylinder doesnt have an end(second) bead and will not be
    /// pushed to cylinder vector, but will be stored in the _pLastCylinder;
    Filament(SubSystem* s, vector<vector<double>>& position, int numBeads, int ID, string projectionType = "STRAIGHT");
    
    ///This destructor is called when a filament is to be removed from the system. Removes all cylinders
    ///and beads associated with the filament
    ~Filament();
    
    void extendFront();  // Addition of a new cylinder (with the first bead = new bead b and last bead = empty: before: --x---x---o, after --x---x---[x---o]). Next position is based on previous beads directions in the filament. This function creates a new bead. So, this function mostly called during further extension, not initiation.
    void extendBack();  // Same as extension front, but adds a new first cylinder with first bead = new bead and a second bead is equal to the firs bead in the cylynder, which used to be first.
    
    ///Extend, used for initialization
    void extendFront(vector<double>& coordinates );
    void extendBack(vector<double>& coordinates );
    
    void retractFront(); // Retraction of front of a cylinder. Removes one cylinder and one bead from the front of filament
    void retractBack(); // Retraction of back of a cylinder. Removes a cylinder and bead from back of filament
    
    ///POLY/DEPOLY
    void polymerizeFront(); //Polymerization of a filament front, which moves the leading bead one monomer length. Updates cylinder parameters accordingly
    void polymerizeBack(); //Same as Polymerization front, but moves the back bead one monomer length, and updates cylinder parameters accordingly
    
    void depolymerizeFront(); // depolymerization of a filament front, which moves the leading bead back one monomer length. Updates cylinder parameters accordingly
    void depolymerizeBack(); // same as depolymerization front, but moves the back bead forward one monomer length. Updates cylinder parameters accordingly
    
    ///Delete a bead from this filament
    void deleteBead(Bead*);
    
    ///More getters
    deque<Cylinder*>& getCylinderVector() {return _cylinderVector;}
    
    ///getters and resetters for deltas
    void resetDeltaPlusEnd() {_deltaPlusEnd = 0;}
    void resetDeltaMinusEnd() {_deltaMinusEnd = 0;}
    
    short getDeltaPlusEnd() {return _deltaPlusEnd;}
    short getDeltaMinusEnd() {return _deltaMinusEnd;}
    
    ///get ID
    int getID() {return _ID;}
    
    //just for example, will rewrite this function, so it not returns anything
    vector<double> nextBeadProjection(Bead* b, double d, vector<double> director);
    vector<vector<double> > straightFilamentProjection(vector<vector<double>>& v, int numBeads);
    vector<vector<double> > zigZagFilamentProjection(vector<vector<double>>& v, int numBeads);
    vector<vector<double> > arcFilamentProjection(vector<vector<double>>& v, int numBeads);
    
    ///Print chemical composition of filament (for debugging only)
    void printChemComposition() {
        for (auto &c : _cylinderVector) {
            c->getCCylinder()->printCCylinder();
        }
    }
};

#endif /* defined(__CytoMech__Filament__) */
