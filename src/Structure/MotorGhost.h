
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v3.0
//
//  Copyright (2015)  Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the MEDYAN web page for more information:
//  http://www.medyan.org
//------------------------------------------------------------------

#ifndef MEDYAN_MotorGhost_h
#define MEDYAN_MotorGhost_h

#include "common.h"

#include "Component.h"
#include "CMotorGhost.h"
#include "MMotorGhost.h"

#include "Database.h"
#include "Histogram.h"
#include "Trackable.h"
#include "Movable.h"
#include "Reactable.h"
#include "RateChangerImpl.h"

//FORWARD DECLARATIONS
class Cylinder;
class Controller;
class DRController;
class SubSystem;

/// A container to store a MMotorGhost and CMotorGhost.
/*!
 *  MotorGhost class is used to manage and store a MMotorGhost and CMotorGhost.
 *  Upon intialization, both of these components are created.
 *
 *  Extending the Movable class, the positions of all instances 
 *  can be updated by the SubSystem.
 *
 *  Extending the Reactable class, the reactions associated with 
 *  all instances can be updated by the SubSystem.
 */
class MotorGhost : public Component, public Trackable, public Movable, public Reactable {
   
friend class Controller;
friend class DRController;

private:
    unique_ptr<MMotorGhost> _mMotorGhost; ///< Pointer to mech motor ghost
    unique_ptr<CMotorGhost> _cMotorGhost; ///< Pointer to chem motor ghost
    
    Cylinder* _c1; ///< First cylinder the motor is bound to
    Cylinder* _c2; ///< Second cylinder the motor is bound to
    
    double _position1; ///< Position on first cylinder
    double _position2; ///< Position on second cylinder
    
    short _motorType; ///< Integer specifying the type of linker
    int _motorID; ///< Integer ID of this motor, managed by Database
    
    float _birthTime; ///< Birth time
    float _walkLength = 0; ///< Walk length of ensemble
    
    Compartment* _compartment; ///< Where this motorghost is
    
    int _numHeads = 1; ///< Number of heads that this motor contains
    double _numBoundHeads = 1; ///< Number of bound heads in the ensemble,
                               ///< which is force-dependent
    
    static Database<MotorGhost*> _motorGhosts;///< Collection in SubSystem
    
    //@{
    ///Histogram data
    static Histogram* _lifetimes;
    static Histogram* _walkLengths;
    //@}
    
    ///For dynamic rate unbinding
    static vector<MotorRateChanger*> _unbindingChangers;
    ///For dynamic rate walking
    static vector<MotorRateChanger*> _walkingChangers;
    
    ///Helper to get coordinate
    void updateCoordinate();

public:
    vector<double> coordinate;
        ///< coordinate of midpoint, updated with updatePosition()
    
    MotorGhost(Cylinder* c1, Cylinder* c2, short motorType,
               double position1 = 0.5, double position2 = 0.5);
    
    virtual ~MotorGhost() noexcept;
    
    //@{
    /// Get cylinder
    Cylinder* getFirstCylinder() {return _c1;}
    Cylinder* getSecondCylinder() {return _c2;}
    //@}
    
    //@{
    /// Set cylinder
    void setFirstCylinder(Cylinder* cylinder) {_c1 = cylinder;}
    void setSecondCylinder(Cylinder* cylinder) {_c2 = cylinder;}
    //@}
    
    /// Set chem motor ghost
    void setCMotorGhost(CMotorGhost* cMotorGhost) {
        _cMotorGhost = unique_ptr<CMotorGhost>(cMotorGhost);
    }
    /// Get chem motor ghost
    CMotorGhost* getCMotorGhost() {return _cMotorGhost.get();}
    
    /// Get mech motor ghost
    MMotorGhost* getMMotorGhost() {return _mMotorGhost.get();}
    
    //@{
    ///Position management function
    double getFirstPosition() {return _position1;}
    void setFirstPosition(double position1) {_position1 = position1;}
    
    double getSecondPosition() {return _position2;}
    void setSecondPosition(double position2) {_position2 = position2;}
    //@}
    
    //@{
    ///Parameter management
    short getType() {return _motorType;}
    int getID() {return _motorID;}
    //@}
    
    /// Get the birth time
    float getBirthTime() {return _birthTime;}
    
    //@{
    /// SubSystem management, inherited from Trackable
    virtual void addToSubSystem() { _motorGhosts.addElement(this);}
    virtual void removeFromSubSystem() {_motorGhosts.removeElement(this);}
    //@}
    
    /// Get all instances of this class from the SubSystem
    static const vector<MotorGhost*>& getMotorGhosts() {
        return _motorGhosts.getElements();
    }
    /// Get the number of motors in this system
    static int numMotorGhosts() {
        return _motorGhosts.countElements();
    }

    /// Get the lifetimes
    static Histogram* getLifetimes() {return _lifetimes;}
    
    /// Get walk lengths
    static Histogram* getWalkLengths() {return _walkLengths;}
    
    /// Update the position, inherited from Movable
    /// @note - changes compartment if needed
    virtual void updatePosition();

    /// Update the reaction rates, inherited from Reactable
    virtual void updateReactionRates();
    
    ///Move a motor head forward
    ///@note - Updates chemical binding and mechanical parameters accordingly
    void moveMotorHead(Cylinder* c, double oldPosition, double newPosition,
                       short boundType, SubSystem* ps);
    
    ///Move a motor head to a new cylinder
    ///@note - Updates chemical binding and mechanical parameters accordingly
    void moveMotorHead(Cylinder* oldC, Cylinder* newC,
                       double oldPosition, double newPosition,
                       short boundType, SubSystem* ps);
    
    virtual void printSelf();
    
    /// Count the number of motor species with a given name in the system
    static species_copy_t countSpecies(const string& name);
};

#endif
