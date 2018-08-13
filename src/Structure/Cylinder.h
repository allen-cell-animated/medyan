
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v3.2
//
//  Copyright (2015-2018)  Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the MEDYAN web page for more information:
//  http://www.medyan.org
//------------------------------------------------------------------

#ifndef MEDYAN_Cylinder_h
#define MEDYAN_Cylinder_h

#include <iostream>

#include "common.h"

#include "MCylinder.h"
#include "CCylinder.h"
#include "RateChanger.h"

#include "Database.h"
#include "Trackable.h"
#include "Movable.h"
#include "Reactable.h"
#include "DynamicNeighbor.h"
#include "Component.h"

//FORWARD DECLARATIONS
class Filament;
class Compartment;

/// A container to store a MCylinder and CCylinder.
/*!
 *  Cylinder class is used to manage and store a MCylinder and CCylinder.
 *  Upon intialization, both of these components are created.
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
class Cylinder : public Component, public Trackable, public Movable,
                                   public Reactable, public DynamicNeighbor {
    
friend class CController;
friend class DRController;
    
private:
    Bead* _b1;  ///< Pointer to the first bead.
    Bead* _b2; ///< Pointer to the end bead.
    
    unique_ptr<MCylinder> _mCylinder; ///< Pointer to mech cylinder
    unique_ptr<CCylinder> _cCylinder; ///< Pointer to chem cylinder
    
    int _position;          ///< Position on structure
    
    bool _plusEnd = false;  ///< If the cylinder is at the plus end
    bool _minusEnd = false; ///< If the cylinder is at the minus end
    
    short _type; ///< Type of cylinder, either corresponding to Filament or other
                                       
    int _ID; ///< Unique ID of cylinder, managed by Database
    
    Compartment* _compartment = nullptr; ///< Where this cylinder is
    
    Cylinder* _branchingCylinder = nullptr; ///< ptr to a branching cylinder
    
    static Database<Cylinder*> _cylinders; ///< Collection in SubSystem
    
    ///For dynamic polymerization rate
    static vector<FilamentRateChanger*> _polyChanger;
                                       
    static ChemManager* _chemManager; ///< A pointer to the ChemManager,
                                      ///< intiailized by CController
    
    ///Helper to get coordinate
    void updateCoordinate();
    
public:
    vector<double> coordinate;
    ///< Coordinates of midpoint, updated with updatePosition()
                                       
    /// Constructor, initializes a cylinder
    Cylinder(Composite* parent, Bead* b1, Bead* b2, short type, int position,
             bool extensionFront = false,
             bool extensionBack  = false,
             bool initialization = false);
                                       
    virtual ~Cylinder() noexcept;
    
    /// Get mech cylinder
    MCylinder* getMCylinder() {return _mCylinder.get();}
    
    /// Get chem cylinder
    CCylinder* getCCylinder() {return _cCylinder.get();}
    /// set chem cylinder
    /// @note: since this is a unique ptr, will implicitly delete old chem cylinder
    void setCCylinder(CCylinder* c) {_cCylinder = unique_ptr<CCylinder>(c);}
    
    /// Get cylinder type
    virtual int getType();
                                       
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
    
    int getPosition() {return _position;}
    
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
                                       
    virtual void printSelf();
                                       
    /// Returns whether a cylinder is within a certain distance from another
    /// Uses the closest point between the two cylinders
    virtual bool within(Cylinder* other, double dist);
};

#endif
