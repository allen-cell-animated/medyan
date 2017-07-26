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
    ///For dynamic polymerization rate // TODO use unique pointers
    array<Cylinder*, 6> _cylinders;
    
    unique_ptr<MCylinder> _mCamkii; ///< Pointer to mech cylinder
    unique_ptr<CCylinder> _cCamkii; ///< Pointer to chem cylinder

    SubSystem* _subSystem; ///< SubSystem pointer

    int _ID; ///< Unique ID of camkii, managed by Database
    
    Compartment* _compartment = nullptr; ///< Where this cylinder is
    
    static Database<Camkii*> _camkiis; ///< Collection in SubSystem
                                          
    static ChemManager* _chemManager; ///< A pointer to the ChemManager,
                                      ///< intiailized by CController
    
    int _position;          ///< Position on structure
    ///Helper to get coordinate
    void updateCoordinate();
    
public:
    vector<double> coordinate;
    ///< Coordinates of midpoint, updated with updatePosition()
                                       
    /// Constructor, initializes a Camkii
    Camkii(Composite* parent, int position, bool initialization = false);
                                       
    virtual ~Camkii();


    /// Get compartment
    Compartment* getCompartment() {return _compartment;}
    
    /// Get ID
    int getID() {return _ID;}

    //@{
    /// SubSystem management, inherited from Trackable
    virtual void addToSubSystem() { _camkiis.addElement(this);}
    virtual void removeFromSubSystem() {_camkiis.removeElement(this);}
    //@}
    
    /// Get all instances of this class from the SubSystem
    static const vector<Camkii*>& getCamkiis() {
        return _camkiis.getElements();
    }

    /// Get the number of Camkii in this system
    static int numCamkiis() {
        return _camkiis.countElements();
    }
    
    /// Update the position, inherited from Movable
    /// @note - changes compartment if needed
    virtual void updatePosition();
    
    /// Update the reaction rates, inherited from Reactable
    virtual void updateReactionRates();

    virtual void printSelf();

    ///Check consistency and correctness of binding sites. Used for debugging.
    virtual bool isConsistent() = 0;
};

#endif
#endif //CAMKII