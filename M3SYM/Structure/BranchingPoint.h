
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

#ifndef M3SYM_BranchingPoint_h
#define M3SYM_BranchingPoint_h

#include "common.h"

#include "MBranchingPoint.h"
#include "CBranchingPoint.h"

#include "Database.h"
#include "Trackable.h"
#include "Movable.h"

//FORWARD DECLARATIONS
class Compartment;
class Cylinder;

/// A container to store a MBranchingPoint and CBranchingPoint.
/*!
 *  BranchingPoint class is used to manage and store a MBranchingPoint and
 *  CBranchingPoint. Upon intialization, both of these components are created.
 *
 *  Extending the Movable class, the positions of all instances 
 *  can be updated by the SubSystem.
 */
class BranchingPoint : public Trackable, public Movable {
    
private:
    unique_ptr<MBranchingPoint> _mBranchingPoint; ///< Pointer to mech branch point
    unique_ptr<CBranchingPoint> _cBranchingPoint; ///< Pointer to chem branch point
    
    Cylinder* _c1; ///< Mother cylinder
    Cylinder* _c2; ///< Branching cylinder
    
    double _position;  ///< Position on mother cylinder
    
    short _branchType; ///< Integer specifying the type
    
    int _branchID;     ///< Integer ID of this specific
                       ///< branch point, managed by the Database
    
    float _birthTime;  ///<Birth time
    
    Compartment* _compartment; ///< Where this branch point is
    
    static Database<BranchingPoint*> _branchingPoints; ///< Collection in SubSystem
    
    ///Helper to get coordinate
    void updateCoordinate();
    
public:
    vector<double> coordinate; ///< coordinate of midpoint,
                               ///< updated with updatePosition()
    
    BranchingPoint(Cylinder* c1, Cylinder* c2,
                   short branchType, double position = 0.5);
    virtual ~BranchingPoint() noexcept;
    
    //@{
    ///Get attached cylinder
    Cylinder* getFirstCylinder() {return _c1;}
    Cylinder* getSecondCylinder() {return _c2;}
    //@}
    
    /// Set chem branch point
    void setCBranchingPoint(CBranchingPoint* cBranchingPoint) {
        _cBranchingPoint = unique_ptr<CBranchingPoint>(cBranchingPoint);
    }
    /// Get chem branch point
    CBranchingPoint* getCBranchingPoint() {return _cBranchingPoint.get();}
    
    /// Get mech branch point
    MBranchingPoint* getMBranchingPoint() {return _mBranchingPoint.get();}
    
    //@{
    /// Position management
    double getPosition() {return _position;}
    void setPosition(double position) {_position = position;}
    //@}
    
    //@{
    /// Get branch parameter
    short getType() {return _branchType;}
    int getID() {return _branchID;}
    //@}
    
    /// Get compartment
    Compartment* getCompartment() {return _compartment;}
    
    /// Get the birth time
    float getBirthTime() {return _birthTime;}
    
    //@{
    /// SubSystem management, inherited from Trackable
    virtual void addToSubSystem() { _branchingPoints.addElement(this);}
    virtual void removeFromSubSystem() {_branchingPoints.removeElement(this);}
    //@}
    
    /// Get all instances of this class from the SubSystem
    static const vector<BranchingPoint*>& getBranchingPoints() {
        return _branchingPoints.getElements();
    }
    /// Get the number of branching points in this system
    static int numBranchingPoints() {
        return _branchingPoints.countElements();
    }

    /// Update the position, inherited from Movable
    /// @note - changes compartment if needed
    virtual void updatePosition();
    
    virtual void printInfo();
    
    /// Count the number of brancher species with a given name in the system
    static species_copy_t countSpecies(const string& name);
};

#endif
