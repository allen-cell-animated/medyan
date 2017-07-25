#ifdef CAMKII
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

#ifndef MEDYAN_Camkii_h
#define MEDYAN_Camkii_h

#include <iostream>
#include "Trackable.h"
#include "Composite.h"
#include "common.h"

#include "Cylinder.h"
#include "Bead.h"
#include "RateChanger.h"

#include "Database.h"
#include "Trackable.h"
#include "Movable.h"
#include "Reactable.h"
#include "DynamicNeighbor.h"
#include "Component.h"

#include <array>

//FORWARD DECLARATIONS
class Filament;
class Compartment;
class Cylinder;
class Bead;

/// A container to store a MCamkii and CCamkii.
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
class Camkii : public Component, public Trackable, public Movable,
                                   public Reactable, public DynamicNeighbor {
    
friend class CController;
friend class DRController; // TODO do we need it?
friend class Controller;

private:
    ///For dynamic polymerization rate
    array<Cylinder*, 6> _cylinders;
    
    unique_ptr<MCylinder> _mCamkii; ///< Pointer to mech cylinder
    unique_ptr<CCylinder> _cCamkii; ///< Pointer to chem cylinder

    SubSystem* _subSystem; ///< SubSystem pointer

    int _ID; ///< Unique ID of cylinder, managed by Database
    
    Compartment* _compartment = nullptr; ///< Where this cylinder is
    
    static Database<Camkii*> _camkiiDB; ///< Collection in SubSystem
                                          
    static ChemManager* _chemManager; ///< A pointer to the ChemManager,
                                      ///< intiailized by CController
    
    int _position;          ///< Position on structure
    ///Helper to get coordinate
    void updateCoordinate();
    
public:
    vector<double> coordinate;
    ///< Coordinates of midpoint, updated with updatePosition()
                                       
    /// Constructor, initializes a Camkii
    Camkii(Composite* parent, int position,
             bool extensionFront = false,
             bool extensionBack  = false,
             bool initialization = false);
                                       
    virtual ~Camkii() noexcept;
    
    /// Get mech cylinder
    MCamkii* getMCamkii() {return _mCamkii.get();}
    
    /// Get chem cylinder
    CCamkii* getCCamkii() {return _cCamkii.get();}
    /// set chem cylinder
    /// @note: since this is a unique ptr, will implicitly delete old chem cylinder
    void setCCylinder(CCylinder* c) {_cCylinder = unique_ptr<CCylinder>(c);}
    
    /// Get cylinder type
    virtual int getType();

    /// Get compartment
    Compartment* getCompartment() {return _compartment;}
    
    /// Get ID
    int getID() {return _ID;}
    
    int getPosition() {return _position;}
    
    // TODO do we need?
    //@{
    /// SubSystem management, inherited from Trackable
    virtual void addToSubSystem() { _camkiiDB.addElement(this);}
    virtual void removeFromSubSystem() {_camkiiDB.removeElement(this);}
    //@}
    
    /// Get all instances of this class from the SubSystem
    static const vector<Camkii*>& getCamkiis() {
        return _camkiiDB.getElements();
    }
    /// Get the number of cylinders in this system
    static int numCamkiis() {
        return _camkiiDB.countElements();
    }
    
    /// Update the position, inherited from Movable
    /// @note - changes compartment if needed
    virtual void updatePosition();
    
    /// Update the reaction rates, inherited from Reactable
    virtual void updateReactionRates();
                                       
    /// Check if this cylinder is grown to full length
    bool isFullLength();
                                       
    virtual void printSelf();
                                       
    /// Returns whether a CaMKII is within a certain distance from 2 cylinders
    /// Uses the closest point between the two cylinders
    virtual bool within(Cylinder* a, Cylinder* b, ,double dist);
};

#endif
#endif //CAMKII