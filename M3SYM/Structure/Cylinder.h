
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

#ifndef M3SYM_Cylinder_h
#define M3SYM_Cylinder_h

#include <iostream>

#include "common.h"

#include "Composite.h"
#include "MCylinder.h"
#include "CCylinder.h"
#include "RateChanger.h"

#include "Database.h"
#include "Trackable.h"
#include "Movable.h"
#include "Reactable.h"
#include "DynamicNeighbor.h"

//FORWARD DECLARATIONS
class Filament;
class Compartment;
class DRController;

/// A container to store a MCylinder and CCylinder.
/*!
 *  Cylinder class is used to manage and store a MCylinder and CCylinder.
 *  Upon intialization, both of these components are created.
 *
 *  Extending the Trackable class, all instances are kept and easily 
 *  accessed by the SubSystem.
 *
 *  Extending the Movable class, the positions of all instances 
 *  can be updated by the SubSystem.
 *
 *  Extending the Reactable class, the reactions associated with 
 *  all instances can be updated by the SubSystem.
 *
 *  Extending the DynamicNeighbor class, all instances can be 
 *  kept in [NeighborLists](@ref NeighborList).
 */
class Cylinder : public Composite, public Trackable, public Movable,
                                   public Reactable, public DynamicNeighbor {
    
friend class DRController;
    
private:
    Bead* _b1;  ///< Pointer to the first bead.
    Bead* _b2; ///< Pointer to the end bead.
    
    unique_ptr<MCylinder> _mCylinder; ///< Pointer to mech cylinder
    unique_ptr<CCylinder> _cCylinder; ///< Pointer to chem cylinder
    
    Filament* _pFilament;  ///< Pointer to parent filament where this cylinder belongs
    int _positionFilament; ///< Position on filament (1st, 2nd, ... etc)
    
    bool _plusEnd = false;  ///< If the cylinder is at the plus end
    bool _minusEnd = false; ///< If the cylinder is at the minus end
    
    int _ID; ///< Unique ID of cylinder, managed by Database
    
    Compartment* _compartment = nullptr; ///< Where this cylinder is
    
    Cylinder* _branchingCylinder = nullptr; ///< ptr to a branching cylinder
    
    static Database<Cylinder*> _cylinders; ///< Collection in SubSystem
    
    ///For dynamic polymerization rate
    static FilamentRateChanger* _polyChanger;
    
    ///Helper to get coordinate
    void updateCoordinate();
    
public:
    vector<double> coordinate;
    ///< Coordinates of midpoint, updated with updatePosition()
                                       
    /// Constructor, initializes a cylinder and
    Cylinder(Filament* f, Bead* b1, Bead* b2, int positionFilament,  
             bool extensionFront = false,
             bool extensionBack = false,
             bool initialization = false);
                                       
    virtual ~Cylinder() noexcept;
    
    /// Get mech cylinder
    MCylinder* getMCylinder() {return _mCylinder.get();}
    
    /// Get chem cylinder
    CCylinder* getCCylinder() {return _cCylinder.get();}
    /// set chem cylinder
    /// @note: since this is a unique ptr, will implicitly delete old chem cylinder
    void setCCylinder(CCylinder* c) {_cCylinder = unique_ptr<CCylinder>(c);}
    
    /// Get parent
    Filament* getFilament() {return _pFilament;}
    /// Set parent
    void setFilament(Filament* pf) {_pFilament = pf;}
    
    //@{
    /// Get beads
    Bead* getFirstBead() {return _b1;}
    Bead* getSecondBead() {return _b2;}
    //@}
    
    //@{
    /// Set beads
    void setFirstBead(Bead* b) {_b1 = b;}
    void setSecondBead(Bead* b) {_b2 = b;}
    //@}
    
    /// Get compartment
    Compartment* getCompartment() {return _compartment;}
    
    //@{
    /// Branching cylinder management
    Cylinder* getBranchingCylinder() {return _branchingCylinder;}
    void setBranchingCylinder(Cylinder* c) {_branchingCylinder = c;}
    //@}
    
    /// Get ID
    int getID() {return _ID;}
    
    ///@{
    /// Set plus and minus end boolean markers
    bool isPlusEnd() {return _plusEnd;}
    void setPlusEnd(bool plusEnd) {_plusEnd = plusEnd;}
    
    bool isMinusEnd() {return _minusEnd;}
    void setMinusEnd(bool minusEnd) {_minusEnd = minusEnd;}
    //@}
    
    int getPositionFilament() {return _positionFilament;}
    
    //@{
    /// SubSystem management, inherited from Trackable
    virtual void addToSubSystem() { _cylinders.addElement(this);}
    virtual void removeFromSubSystem() {_cylinders.removeElement(this);}
    //@}
    
    /// Get all instances of this class from the SubSystem
    static const vector<Cylinder*>& getCylinders() {
        return _cylinders.getElements();
    }
    /// Get the number of cylinders in this system
    static int numCylinders() {
        return _cylinders.countElements();
    }
    
    /// Update the position, inherited from Movable
    /// @note - changes compartment if needed
    virtual void updatePosition();
    
    /// Update the reaction rates, inherited from Reactable
    virtual void updateReactionRates();
                                       
    /// Check if this cylinder is grown to full length
    bool isFullLength();

};

#endif
