
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v4.0
//
//  Copyright (2015-2018)  Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the MEDYAN web page for more information:
//  http://www.medyan.org
//------------------------------------------------------------------

#ifndef MEDYAN_BranchingPoint_h
#define MEDYAN_BranchingPoint_h

#include "common.h"

#include "MBranchingPoint.h"
#include "CBranchingPoint.h"

#include "Database.h"
#include "Trackable.h"
#include "Movable.h"
#include "Component.h"
#include "Reactable.h"
#include "RateChangerImpl.h"

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
class BranchingPoint : public Component, public Trackable, public Movable, public Reactable {
    
    friend class Controller;
    friend class DRController;
    
private:
    unique_ptr<MBranchingPoint> _mBranchingPoint; ///< Pointer to mech branch point
    unique_ptr<CBranchingPoint> _cBranchingPoint; ///< Pointer to chem branch point
    
    Cylinder* _c1; ///< Mother cylinder
    Cylinder* _c2; ///< Branching cylinder
    
    floatingpoint _position;  ///< Position on mother cylinder
    
    short _branchType; ///< Integer specifying the type
    
    int _branchID;     ///< Integer ID of this specific
                       ///< branch point, managed by the Database
    
    float _birthTime;  ///<Birth time
    
    Compartment* _compartment; ///< Where this branch point is
    
    static Database<BranchingPoint*> _branchingPoints; ///< Collection in SubSystem
    
    ///Helper to get coordinate
    void updateCoordinate();

    ///For dynamic rate unbinding
    static vector<BranchRateChanger*> _unbindingChangers;
    
public:
    vector<floatingpoint> coordinate; ///< coordinate of midpoint,
                               ///< updated with updatePosition()
    
    BranchingPoint(Cylinder* c1, Cylinder* c2,
                   short branchType, floatingpoint position = 0.5);
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
    floatingpoint getPosition() {return _position;}
    void setPosition(floatingpoint position) {_position = position;}
    //@}
    
    //@{
    /// Get branch parameter
    virtual int getType() {return _branchType;}
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
    
    virtual void printSelf();
    
    /// Count the number of brancher species with a given name in the system
    static species_copy_t countSpecies(const string& name);
    
    /// Update the position, inherited from Movable
    /// @note - changes compartment if needed
    virtual void updatePosition();
    
    //Qin ------
    /// Update the reaction rates, inherited from Reactable
    virtual void updateReactionRates();


};

#endif
