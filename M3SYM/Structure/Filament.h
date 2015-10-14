
//------------------------------------------------------------------
//  **M3SYM** - Simulation Package for the Mechanochemical
//              Dynamics of Active Networks, 3rd Generation
//
//  Copyright (2014) Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the Papoian lab page for more information:
//  http://papoian.chem.umd.edu/
//------------------------------------------------------------------

#ifndef M3SYM_Filament_h
#define M3SYM_Filament_h

#include <vector>
#include <iostream>
#include <deque>
#include <cmath>

#include "common.h"

#include "Database.h"
#include "Trackable.h"

//FORWARD DECLARATIONS
class SubSystem;
class Cylinder;
class Bead;

/// Used to store data about connectivity of [Cylinders](@ref Cylinder) and [Beads] (@ref Bead).
/*!
 * This class contains information about [Cylinders](@ref Cylinder) and [Beads]
 * (@ref Bead) connectivity. It has functionality to polymerize and depolymerize,
 * as well as extend and retract by creating and deleting Cylinder and Bead in the 
 * system.
 *
 * Extending the Trackable class, all instances are kept and 
 * easily accessed by the SubSystem.
 * 
 * A Filament can also be initialized as a number of different shapes.
 */
class Filament : public Trackable {

private:
    short _filType; ///< Filament type
    
    deque<Cylinder*> _cylinderVector; ///< Vector of cylinders;
    SubSystem* _subSystem; ///< SubSystem pointer
    
    int _ID; ///< Unique integer id of this filament
    
    short _deltaPlusEnd = 0;   ///< Change in filament's cylinders
                               ///< at plus end since last snapshot
    short _deltaMinusEnd = 0;  ///< Change in filament's cylinders
                               ///< at minus end since last snapshot
    
    static Database<Filament*> _filaments; ///< Collection in SubSystem
    
public:
    /// This constructor creates a short filament, containing only two beads, at runtime.
    /// Coordinates of the first bead is an input, second is set up by using an input
    /// direction. Using all this, two constructors for beads and cylinders are called.
    /// @param nucleation - this filament was nucleated at runtime by a non-branching species
    /// @param branching - this filament was branched at runtime from an existing filament
	Filament(SubSystem* s, short filamentType,
                           vector<double>& position,
                           vector<double>& direction,
                           bool nucleation = false,
                           bool branch = false);
    
    /// This constructor is called to create a filament at startup. It creates a filament
    /// with a number of beads numBeads. Filaments starts and ends in the point
    /// determined by position vector.
    Filament(SubSystem* s, short filamentType,
             vector<vector<double>>& position, int numBeads,
             string projectionType = "STRAIGHT");
    
    /// This constructor is called when a filament is severed. It creates a filament
    /// that initially has no cylinders.
    Filament(SubSystem* s, short filamentType)
        : Trackable(), _subSystem(s), _filType(filamentType), _ID(_filaments.getID()) {}
    
    /// This destructor is called when a filament is to be removed from the system.
    /// Removes all cylinders and beads associated with the filament.
    ~Filament();
    
    /// Addition of a new cylinder. Next position is based on previous beads directions
    /// in the filament. This function creates a new bead. So, this function is mostly
    /// called during further extension, not initiation.
    /// @param plusEnd - the plus end species to be marked. Only if doing chemistry.
    void extendFront(short plusEnd);
    /// Same as extension front, but adds a new first cylinder with first bead = new
    /// bead and a second bead is equal to the first bead in the cylinder, which used
    /// to be first.
    /// @param minusEnd - the minus end species to be marked. Only if doing chemistry.
    void extendBack(short minusEnd);
    
    ///Extend, used for initialization
    void extendFront(vector<double>& coordinates);
    void extendBack(vector<double>& coordinates);
    
    /// Retraction of front of a cylinder. Removes one cylinder and one bead from the
    /// front of filament.
    void retractFront();
    /// Retraction of back of a cylinder. Removes a cylinder and bead from back of
    /// filament.
    void retractBack();
    
    /// Polymerization of a filament front, which moves the leading bead one monomer
    /// length. Updates cylinder parameters accordingly.
    void polymerizeFront();
    /// Same as Polymerization front, but moves the back bead one monomer length, and
    /// updates cylinder parameters accordingly.
    void polymerizeBack();
    
    /// Depolymerization of a filament front, which moves the leading bead back one
    /// monomer length. Updates cylinder parameters accordingly.
    void depolymerizeFront();
    /// Same as depolymerization front, but moves the back bead forward one monomer
    /// length. Updates cylinder parameters accordingly.
    void depolymerizeBack();
    
    /// Initialize the nucleation of a new filament
    /// Initializes all chemical species in initial Cylinder
    void nucleate(short plusEnd, short filament, short minusEnd);
    
    /// Sever a filament at a given cylinder position. The back part of the filament
    /// will remain as this (original) filament, and the new filament, which is
    /// the front part of the originally severed filament, will be returned.
    Filament* sever(int cylinderPosition);
    
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
    
    /// Get type
    short getType() {return _filType;}
    
    //@{
    /// SubSystem management, inherited from Trackable
    virtual void addToSubSystem() { _filaments.addElement(this);}
    virtual void removeFromSubSystem() {_filaments.removeElement(this);}
    //@}
    
    /// Get all instances of this class from the SubSystem
    static const vector<Filament*>& getFilaments() {
        return _filaments.getElements();
    }
    /// Get the number of filaments in this system
    static int numFilaments() {
        return _filaments.countElements();
    }
    
    //@{
    /// Projection function, returns a vector of coordinates for bead creation
    vector<double> nextBeadProjection(Bead* b, double d, vector<double> director);
    
    vector<vector<double>>
    straightFilamentProjection(vector<vector<double>>& v, int numBeads);
    
    vector<vector<double>>
    zigZagFilamentProjection(vector<vector<double>>& v, int numBeads);
    
    vector<vector<double>>
    arcFilamentProjection(vector<vector<double>>& v, int numBeads);
    //@}
    
    virtual void printInfo();
    
    /// Check the consistency of the filament. For now,
    /// Mainly involves checking chemistry of the constituent cylinders.
    bool isConsistent();
    
    /// Count the number of filament species with a given name in the system
    static species_copy_t countSpecies(short filamentType, const string& name);
};

#endif
