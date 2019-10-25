
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
class MotorGhost : public Component, public Trackable, public Movable, public Reactable,
    public Database< MotorGhost, false > {
   
friend class Controller;
friend class DRController;
    
friend class MotorBindingManager;
    
friend struct UpdateMotorIDCallback;
friend struct MotorBindingCallback;
friend struct MotorUnbindingCallback;

#ifdef COLLECTMOTORDATA
struct MotorGhostData{
    vector<int> ID;
    vector<int> Nwalkingsteps;
    vector<float> Lifetime;
    vector<float> Energyvec;
};

static vector<MotorGhostData> _motorGhostdatavec;
#endif
private:

    chrono::high_resolution_clock::time_point mins, mine;
    unique_ptr<MMotorGhost> _mMotorGhost; ///< Pointer to mech motor ghost
    unique_ptr<CMotorGhost> _cMotorGhost; ///< Pointer to chem motor ghost
    
    Cylinder* _c1; ///< First cylinder the motor is bound to
    Cylinder* _c2; ///< Second cylinder the motor is bound to
    
    floatingpoint _position1; ///< Position on first cylinder
    floatingpoint _position2; ///< Position on second cylinder
    
    short _motorType; ///< Integer specifying the type of linker
    
    float _birthTime; ///< Birth time
    float _walkLength = 0; ///< Walk length of ensemble
    
    Compartment* _compartment; ///< Where this motorghost is
    
    int _numHeads = 1; ///< Number of heads that this motor contains
    floatingpoint _numBoundHeads = 1; ///< Number of bound heads in the ensemble,
                               ///< which is force-dependent
    
    //@{
    ///Kinetic rates of individual motor heads
    floatingpoint _onRate = 0.0;
    floatingpoint _offRate = 0.0;
    //@}
    
    //@{
    ///Histogram data
    static Histogram* _lifetimes;
    static Histogram* _walkLengths;
    //@}
    
    ///For dynamic rate unbinding
    static vector<MotorRateChanger*> _unbindingChangers;
    ///For dynamic rate walking
    static vector<MotorRateChanger*> _walkingChangers;



public:
    #ifdef MOTORBIASCHECK
    size_t walkingsteps = 0;
    #endif
    vector<floatingpoint> coordinate;
        ///< coordinate of midpoint, updated with updatePosition()
    
    ///Standard constructor
    MotorGhost(Cylinder* c1, Cylinder* c2, short motorType,
               floatingpoint position1 = 0.5, floatingpoint position2 = 0.5,
               floatingpoint onRate = 0.0, floatingpoint offRate = 0.0);
    
    virtual ~MotorGhost() noexcept;
    
    ///Helper to get coordinate
    void updateCoordinate();
    
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
    floatingpoint getFirstPosition() {return _position1;}
    void setFirstPosition(floatingpoint position1) {_position1 = position1;}
    
    floatingpoint getSecondPosition() {return _position2;}
    void setSecondPosition(floatingpoint position2) {_position2 = position2;}
    //@}
    
    //@{
    ///Parameter management
    virtual int getType() {return _motorType;}
    //@}
    
    /// Get the birth time
    float getBirthTime() {return _birthTime;}
    
    // get num heads
    int getNumHeads(){return _numHeads;}
    
    //@{
    /// SubSystem management, inherited from Trackable
    // Does nothing
    virtual void addToSubSystem() override { }
    virtual void removeFromSubSystem() override {}
    //@}
    
    /// Get all instances of this class from the SubSystem
    static const vector<MotorGhost*>& getMotorGhosts() {
        return getElements();
    }
    /// Get the number of motors in this system
    static int numMotorGhosts() {
        return getElements().size();
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
    void moveMotorHead(Cylinder* c, floatingpoint oldPosition, floatingpoint newPosition,
                       short boundType, SubSystem* ps);
    
    ///Move a motor head to a new cylinder
    ///@note - Updates chemical binding and mechanical parameters accordingly
    void moveMotorHead(Cylinder* oldC, Cylinder* newC,
                       floatingpoint oldPosition, floatingpoint newPosition,
                       short boundType, SubSystem* ps);
    
    virtual void printSelf();
    
    /// Count the number of motor species with a given name in the system
    static species_copy_t countSpecies(const string& name);
};

#endif
