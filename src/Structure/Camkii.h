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

#include "Database.h"
#include "Trackable.h"
#include "Movable.h"
#include "Reactable.h"
#include "DynamicNeighbor.h"

#include "ChemManager.h"

#include <array>

//FORWARD DECLARATIONS
class Filament;
class Compartment;
class Cylinder;
class Bead;
class SubSystem;

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
class Camkii : public Composite, public Trackable, public Movable,
                                   public Reactable, public DynamicNeighbor {
    
friend class CController;
friend class DRController; // TODO do we need it?
friend class Controller;

private:
    ///For dynamic polymerization rate // TODO use unique pointers
    array<Cylinder*, 6> _cylinders;

    short _type; ///< Filament type

    SubSystem* _subSystem; ///< SubSystem pointer

    int _ID; ///< Unique ID of camkii, managed by Database
    
    Compartment* _compartment = nullptr; ///< Where this cylinder is
    
    static Database<Camkii*> _camkiiDB; ///< Collection in SubSystem
                                          
    static ChemManager* _chemManager; ///< A pointer to the ChemManager,
                                      ///< intiailized by CController

    void hexagonProjection(vector<double>& filamentInterfacePoint);
    ///Helper to get coordinate
    void updateCoordinate();
    
public:
    vector<double> coordinate;
    ///< Coordinates of midpoint, updated with updatePosition()
                                       
    /// Constructor, initializes a Camkii
    Camkii(SubSystem* s, short type, vector<double>& filamentInterfacePoint);
                                       
    virtual ~Camkii();


    /// Get compartment
    Compartment* getCompartment() {return _compartment;}
    
    /// Get ID
    int getID() {return _ID;}

    //@{
    /// SubSystem management, inherited from Trackable
    virtual void addToSubSystem() { _camkiiDB.addElement(this);}
    virtual void removeFromSubSystem() {_camkiiDB.removeElement(this);}
    //@}
    
    /// Get all instances of this class from the SubSystem
    static const vector<Camkii*>& getCamkiis() {
        return _camkiiDB.getElements();
    }

    /// Get the number of Camkii in this system
    static int numCamkiis() {
        return _camkiiDB.countElements();
    }
    
    /// Update the position, inherited from Movable
    /// @note - changes compartment if needed
    virtual void updatePosition();
    
    /// Update the reaction rates, inherited from Reactable
    virtual void updateReactionRates();

    virtual void printSelf();

    ///Check consistency and correctness of binding sites. Used for debugging.
    virtual bool isConsistent();
                                       
    virtual int getType() {return _type;}
};

#endif
#endif //CAMKII
