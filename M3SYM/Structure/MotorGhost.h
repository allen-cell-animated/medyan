
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

#ifndef M3SYM_MotorGhost_h
#define M3SYM_MotorGhost_h

#include "common.h"

#include "MotorGhostDB.h"

#include "Composite.h"
#include "CMotorGhost.h"
#include "MMotorGhost.h"
#include "Movable.h"
#include "Reactable.h"
#include "RateChangerImpl.h"

#include "SysParams.h"
//FORWARD DECLARATIONS
class Cylinder;
class DRController;

/// A container to store a MMotorGhost and CMotorGhost.
/*!
 * MotorGhost class is used to manage and store a MMotorGhost and CMotorGhost.
 * Upon intialization, both of these components are created. Extending the Movable and
 * Reactable classes, the MotorGhost can update its position and reactions according to 
 * mechanical equilibration.
 */
class MotorGhost : public Composite, public Movable, public Reactable {
   
friend class DRController;

private:
    unique_ptr<MMotorGhost> _mMotorGhost; ///< Pointer to mech motor ghost
    unique_ptr<CMotorGhost> _cMotorGhost; ///< Pointer to chem motor ghost
    
    Cylinder* _c1; ///< First cylinder the linker is bound to
    Cylinder* _c2; ///< Second cylinder the linker is bound to
    
    double _position1; ///< Position on first cylinder
    double _position2; ///< Position on second cylinder
    
    short _motorType; ///< Integer specifying the type of linker
    int _motorID; ///< Integer ID of this motor, managed by MotorGhostDB
    
    float _birthTime; ///< Birth time
    
    Compartment* _compartment; ///< Where this motorghost is
    
    int _numHeads = 1; ///< Number of heads that this motor contains
    
    ///For dynamic rate unbinding
    static vector<MotorRateChanger*> _unbindingChangers;
    ///For dynamic rate walking
    static vector<MotorRateChanger*> _walkingChangers;

public:
    vector<double> coordinate;
        ///< coordinate of midpoint, updated with updatePosition()
    
    MotorGhost(Cylinder* c1, Cylinder* c2, short motorType,
               double position1 = 0.5, double position2 = 0.5, bool creation = false);
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
    short getMotorType() {return _motorType;}
    int getMotorID() {return _motorID;}
    //@}
    
    /// Get the birth time
    float getBirthTime() {return _birthTime;}
    
    /// Update the position
    /// @note - changes compartment if needed
    virtual void updatePosition();

    /// Update the reaction rates
    virtual void updateReactionRates();

};

#endif
