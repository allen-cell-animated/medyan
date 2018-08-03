
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v3.1
//
//  Copyright (2015-2016)  Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the MEDYAN web page for more information:
//  http://www.medyan.org
//------------------------------------------------------------------

#ifndef MEDYAN_Filament_h
#define MEDYAN_Filament_h

#include <vector>
#include <iostream>
#include <deque>
#include <cmath>

#include "common.h"

#include "Database.h"
#include "Histogram.h"
#include "Trackable.h"
#include "Composite.h"

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
class Filament : public Composite, public Trackable {

friend class Controller;
    
private:
    /// Deque of cylinders
    /// @note - the "front" of this deck is the minus end of the filament.
    deque<Cylinder*> _cylinderVector;
    
    SubSystem* _subSystem; ///< SubSystem pointer
    
    short _filType; ///< Filament type
    
    int _ID; ///< Unique integer id of this filament
    
    short _deltaPlusEnd = 0;   ///< Change in filament's cylinders
                               ///< at plus end since last snapshot
    short _deltaMinusEnd = 0;  ///< Change in filament's cylinders
                               ///< at minus end since last snapshot
    
    int _plusEndPosition   = 0;  ///< Position of plus end bead at last turnover
    double _turnoverTime   = 0;  ///< Time since last turnover
    
    static Database<Filament*> _filaments; ///< Collection in SubSystem
    
    //@{
    ///Histogram data
    static Histogram* _turnoverTimes;
    //@}
    
public:
    /// This constructor creates a short filament, containing only two beads, at runtime.
    /// Coordinates of the first bead is an input, second is set up by using an input
    /// direction. Using all this, two constructors for beads and cylinders are called.
    /// @param nucleation - this filament was nucleated at runtime by a non-camkiiing species
    /// @param camkiiing - this filament was camkiied at runtime from an existing filament
	Filament(SubSystem* s, short filamentType,
                           vector<double>& position,
                           vector<double>& direction,
                           bool nucleation = false,
                           bool camkii = false);
    
    /// This constructor is called to create a filament at startup. It creates a filament
    /// with a number of beads numBeads. Filaments starts and ends in the point
    /// determined by position vector.
    Filament(SubSystem* s, short filamentType,
             vector<vector<double>>& position, int numBeads,
             string projectionType = "PREDEFINED");
    
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
    void extendPlusEnd(short plusEnd);
    /// Same as extension front, but adds a new first cylinder with first bead = new
    /// bead and a second bead is equal to the first bead in the cylinder, which used
    /// to be first.
    /// @param minusEnd - the minus end species to be marked. Only if doing chemistry.
    void extendMinusEnd(short minusEnd);
    
    ///Extend, used for initialization
    void extendPlusEnd(vector<double>& coordinates);
    void extendMinusEnd(vector<double>& coordinates);
    
    /// Retraction of front of a cylinder. Removes one cylinder and one bead from the
    /// front of filament.
    void retractPlusEnd();
    /// Retraction of back of a cylinder. Removes a cylinder and bead from back of
    /// filament.
    void retractMinusEnd();
    
    /// Polymerization of a filament front, which moves the leading bead one monomer
    /// length. Updates cylinder parameters accordingly.
    void polymerizePlusEnd();
    /// Same as Polymerization front, but moves the back bead one monomer length, and
    /// updates cylinder parameters accordingly.
    void polymerizeMinusEnd();
    
    /// Depolymerization of a filament front, which moves the leading bead back one
    /// monomer length. Updates cylinder parameters accordingly.
    void depolymerizePlusEnd();
    /// Same as depolymerization front, but moves the back bead forward one monomer
    /// length. Updates cylinder parameters accordingly.
    void depolymerizeMinusEnd();
    
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
    int getType() {return _filType;}
    
    //@{
    /// Get end cylinder
    Cylinder* getMinusEndCylinder() {return _cylinderVector.front();}
    Cylinder* getPlusEndCylinder() {return _cylinderVector.back();}
    //@}
    
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
    /// Get the turnover times
    static Histogram* getTurnoverTimes() {return _turnoverTimes;}
    
    //@{
    /// Projection function, returns a vector of coordinates for bead creation
    vector<double> nextBeadProjection(Bead* b, double d, vector<double> director);
    
    vector<vector<double>> straightFilamentProjection(vector<vector<double>>& v, int numBeads);
    vector<vector<double>> zigZagFilamentProjection(vector<vector<double>>& v, int numBeads);
    vector<vector<double>> arcFilamentProjection(vector<vector<double>>& v, int numBeads);
    //Aravind 18 Feb 2016.
    vector<vector<double>> predefinedFilamentProjection(vector<vector<double>>& v, int numBeads);
    //@}
    
    virtual void printSelf();
    
    /// Check the consistency of the filament. For now,
    /// Mainly involves checking chemistry of the constituent cylinders.
    bool isConsistent();
    
    /// Count the number of filament species with a given name in the system
    static species_copy_t countSpecies(short filamentType, const string& name);
};

#endif