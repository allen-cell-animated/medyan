
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

#ifndef MEDYAN_CaMKIIingPoint_h
#define MEDYAN_CaMKIIingPoint_h

#include "common.h"

#include "MCaMKIIingPoint.h"
#include "CCaMKIIingPoint.h"

#include "Database.h"
#include "Trackable.h"
#include "Movable.h"
#include "Component.h"
#include "CaMKIICylinder.h"

//FORWARD DECLARATIONS
class Compartment;
class Cylinder;

/// A container to store a MCaMKIIingPoint and CCaMKIIingPoint.
/*!
 *  CaMKIIingPoint class is used to manage and store a MCaMKIIingPoint and
 *  CCaMKIIingPoint. Upon intialization, both of these components are created.
 *
 *  Extending the Movable class, the positions of all instances 
 *  can be updated by the SubSystem.
 */
class CaMKIIingPoint : public Component, public Trackable, public Movable {
    
private:
    unique_ptr<MCaMKIIingPoint> _mCaMKIIingPoint; ///< Pointer to mech camkii point
	  unique_ptr<CCaMKIIingPoint> _cCaMKIIingPoint; ///< Pointer to chem camkii point
    unique_ptr<CaMKIICylinder> _camkiiCylinder; ///< Pointer to CaMKIICylinder


    vector<tuple<Cylinder*, short>> _bonds;

    short _camkiiType; ///< Integer specifying the type
    
    int _camkiiID;     ///< Integer ID of this specific
                       ///< camkii point, managed by the Database
    
    float _birthTime;  ///<Birth time
    const short _filType = CAMKII_CYLINDER_FILAMENT_TYPE;

    Compartment* _compartment; ///< Where this camkii point is
    
    static Database<CaMKIIingPoint*> _camkiiingPoints; ///< Collection in SubSystem

    ///Helper to get coordinate
    void updateCoordinate(vector<double> *coordinate);
    
    void updateCaMKIIingPointCoM();


public:
	/*
	 * The "coordinate" in CaMKIIingPoint is a pointer
	 * while
	 * the "coordinate" in BranchingPoint is a value.
	 *
	 * This coordinate basically points to the bead that creates the CaMKIIingPoint.
	 */
    vector<double> *coordinate; ///< coordinate of midpoint,
                               ///< updated with updatePosition()
    
    CaMKIIingPoint(Cylinder* cylinder, short camkiiType, double position);
    virtual ~CaMKIIingPoint()  noexcept;

	//@{
	///Get attached bonds tuples <cylinders, short> (a cylinder and a position)
	tuple<Cylinder *, short> getBond(int n) { return _bonds.at(n); }

	Cylinder *getCylinder(int n) { return get<0>(getBond(n)); }

	void addBond(Cylinder *c, short pos);

	Cylinder *getFirstCylinder() { return get<0>(_bonds.at(0)); }

	Cylinder *getSecondCylinder() { return get<0>(_bonds.at(1)); } //TODO fix
	int getCoordinationNumber() { return _bonds.size(); }

	CaMKIICylinder *getCaMKIICylinder() { return _camkiiCylinder.get(); }

	void setCaMKIICylinder(CaMKIICylinder *CaMKIIcylinder) {
		_camkiiCylinder = unique_ptr<CaMKIICylinder>(CaMKIIcylinder);
	}

	vector<tuple<Cylinder *, short>> &getBonds() { return _bonds; }
	//@}
    
    
    /// Set chem camkii point
    void setCCaMKIIingPoint(CCaMKIIingPoint* cCaMKIIingPoint) {
        _cCaMKIIingPoint = unique_ptr<CCaMKIIingPoint>(cCaMKIIingPoint);
    }
    /// Get chem camkii point
    CCaMKIIingPoint* getCCaMKIIingPoint() {return _cCaMKIIingPoint.get();}
    
    /// Get mech camkii point
    MCaMKIIingPoint* getMCaMKIIingPoint() {return _mCaMKIIingPoint.get();}

    //@{
    /// Position management
    short getPosition() {return get<1>(_bonds.at(0));} //TODO fix later
    //void setPosition(double position) {_position = position;}
    //@}
    
    //@{
    /// Get camkii parameter
    virtual int getType() {return _camkiiType;}
    int getID() {return _camkiiID;}
    //@}
    
    /// Get compartment
    Compartment* getCompartment() {return _compartment;}
    
    /// Get the birth time
    float getBirthTime() {return _birthTime;}
    
    //@{
    /// SubSystem management, inherited from Trackable
    virtual void addToSubSystem() { _camkiiingPoints.addElement(this);}
    virtual void removeFromSubSystem() {_camkiiingPoints.removeElement(this);}
    //@}

    /// Get all instances of this class from the SubSystem
    static const vector<CaMKIIingPoint*>& getCaMKIIingPoints() {
        return _camkiiingPoints.getElements();
    }
    /// Get the number of camkiiing points in this system
    static int numCaMKIIingPoints() {
        return _camkiiingPoints.countElements();
    }

    /// Update the position, inherited from Movable
    /// @note - changes compartment if needed
    virtual void updatePosition();
    
    virtual void printSelf();
    
    /// Count the number of camkiier species and dummy species with a given name in the system
	static species_copy_t countSpecies(const string& name);
	static species_copy_t countDummySpecies(const string& name);

	void removeRandomBond();

    void updateReactionRates();
};

#endif
