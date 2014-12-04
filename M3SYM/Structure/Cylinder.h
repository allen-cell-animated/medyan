
//------------------------------------------------------------------
//  **M3SYM** - Simulation Package for the Mechanochemical
//              Dynamics of Active Networks, 3rd Generation
//
//  Copyright (2014) Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the Papoian lab page for installation and documentation:
//  http://papoian.chem.umd.edu/
//------------------------------------------------------------------

#ifndef M3SYM_Cylinder_h
#define M3SYM_Cylinder_h

#include <iostream>

#include "common.h"

#include "CylinderDB.h"

#include "Composite.h"
#include "MCylinder.h"
#include "CCylinder.h"
#include "Neighbor.h"

#include "Movable.h"
#include "Reactable.h"

//FORWARD DECLARATIONS
class Filament;
class Compartment;

/// A container to store a MCylinder and CCylinder.
/*!
 * Cylinder class is used to manage and store a MCylinder and CCylinder.
 * Upon intialization, both of these components are created. Extending the Movable and Reactable
 * classes, the Cylinder can update its position and reactions according to mechanical equilibration.
 */

class Cylinder : public Composite, public Neighbor, public Movable, public Reactable {
    
private:
    Bead* _b1;  ///< Pointer to the first bead.
    Bead* _b2; ///< Pointer to the end bead.
    
    unique_ptr<MCylinder> _mCylinder; ///< Pointer to mech cylinder
    unique_ptr<CCylinder> _cCylinder; ///< Pointer to chem cylinder
    
    Filament* _pFilament; //< Pointer to parent filament where this cylinder belongs;
    int _positionFilament; ///< Position on filament (1st, 2nd, ... etc)
    bool _last = false; ///< If the cylinder is last in the filament's cylinder list
    
    int _ID; ///< Unique ID of cylinder, managed by CylinderDB
    
    Compartment* _compartment = nullptr; ///< Where this cylinder is
    
public:
    vector<double> coordinate; ///< Coordinates of midpoint, updated with updatePosition()
    
    Cylinder(Filament* f, Bead* b1, Bead* b2, int positionFilament,  
             bool extensionFront = false, bool extensionBack = false, bool creation = false);
    ~Cylinder();
    
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
    
    /// Get compartment
    Compartment* getCompartment() {return _compartment;}
    
    /// Get ID
    const int getID() {return _ID;}
    
    bool last(){ return _last;}
    void setLast(bool b){ _last = b;}
    
    int getPositionFilament() {return _positionFilament;}
    
    /// Update the position
    /// @note - changes compartment if needed
    virtual void updatePosition();
    
    /// Update the reaction rates
    virtual void updateReactionRates();

};

#endif
