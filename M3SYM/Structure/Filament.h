
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

#ifndef M3SYM_Filament_h
#define M3SYM_Filament_h

#include <vector>
#include <iostream>
#include <deque>
#include <cmath>

#include "common.h"

#include "FilamentDB.h"

//FORWARD DECLARATIONS
class SubSystem;
class Cylinder;
class Bead;

/// Used to store data about connectivity of [Cylinders](@ref Cylinder) and [Beads](@ref Bead).
/*!
 * This class contains information about [Cylinders](@ref Cylinder) and [Beads](@ref Bead) connectivity. It has
 * functionality to polymerize and depolymerize, as well as extend and retract by creating and deleting Cylinders and 
 * Beads in the system. 
 * 
 * A Filament can also be initialized as a number of shapes, including a zig zag and arc projection.
 */
class Filament {

private:
    deque<Cylinder*> _cylinderVector; ///< Vector of cylinders;
    SubSystem* _subSystem; ///< SubSystem pointer
    
    int _ID; ///< Unique integer id of this filament
    
    short _deltaPlusEnd = 0;   ///< Change in filament's cylinders at plus end since last snapshot
    short _deltaMinusEnd = 0; ///< Change in filament's cylinders at minus end since last snapshot
    
public:
    /// This constructor creates a short filament, containing only two beads. Coordinates of the first bead is an
    /// input, second is set up by using an input direction. Using all this, two constructors
    /// for beads and cylinders are called.
	Filament(SubSystem* s, vector<double>& position, vector<double>& direction);
    
    /// This constructor is called to create a longer filament. It creates a filament with a number of beads numBeads.
    /// Filaments starts and ends in the point determined by position vector and has a direction direction. Number of
    /// beads is equal to the number of cylinders. The last cylinder doesnt have an end(second) bead and will not be
    /// pushed to cylinder vector, but will be stored in the _pLastCylinder;
    Filament(SubSystem* s, vector<vector<double>>& position, int numBeads, string projectionType = "STRAIGHT");
    
    ///This destructor is called when a filament is to be removed from the system. Removes all cylinders
    ///and beads associated with the filament.
    ~Filament();
    
    /// Addition of a new cylinder (with the first bead = new bead b and last bead = empty: before: --x---x---o, after --x---x---[x---o]).
    /// Next position is based on previous beads directions in the filament. This function creates a new bead. So, this function mostly called
    /// during further extension, not initiation.
    void extendFront();
    /// Same as extension front, but adds a new first cylinder with first bead = new bead and a second bead is equal to the firs bead in the
    /// cylynder, which used to be first.
    void extendBack();
    
    ///Extend, used for initialization
    void extendFront(vector<double>& coordinates );
    void extendBack(vector<double>& coordinates );
    
    /// Retraction of front of a cylinder. Removes one cylinder and one bead from the front of filament.
    void retractFront();
    /// Retraction of back of a cylinder. Removes a cylinder and bead from back of filament.
    void retractBack();
    
    /// Polymerization of a filament front, which moves the leading bead one monomer length. Updates cylinder parameters accordingly.
    void polymerizeFront();
    /// Same as Polymerization front, but moves the back bead one monomer length, and updates cylinder parameters accordingly.
    void polymerizeBack();
    
    /// Depolymerization of a filament front, which moves the leading bead back one monomer length. Updates cylinder parameters accordingly.
    void depolymerizeFront();
    /// Same as depolymerization front, but moves the back bead forward one monomer length. Updates cylinder parameters accordingly.
    void depolymerizeBack();
    
    /// Get vector of cylinders that this filament contains.
    deque<Cylinder*>& getCylinderVector() {return _cylinderVector;}
    
    //@{
    /// Reset delta
    void resetDeltaPlusEnd() {_deltaPlusEnd = 0;}
    void resetDeltaMinusEnd() {_deltaMinusEnd = 0;}
    //@}
    
    //@{
    /// Get delta
    short getDeltaPlusEnd() {return _deltaPlusEnd;}
    short getDeltaMinusEnd() {return _deltaMinusEnd;}
    //@}
    
    /// Get ID
    int getID() {return _ID;}
    
    //@{
    /// Projection function, returns a vector of coordinates for bead creation
    vector<double> nextBeadProjection(Bead* b, double d, vector<double> director);
    vector<vector<double> > straightFilamentProjection(vector<vector<double>>& v, int numBeads);
    vector<vector<double> > zigZagFilamentProjection(vector<vector<double>>& v, int numBeads);
    vector<vector<double> > arcFilamentProjection(vector<vector<double>>& v, int numBeads);
    //@}
    
    /// Print chemical composition of filament (for debugging only)
    void printChemComposition();
};

#endif
