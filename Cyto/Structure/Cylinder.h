//
//  Cylinder.h
//  Cyto
//
//  Created by James Komianos on 7/31/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#ifndef __Cyto__Cylinder__
#define __Cyto__Cylinder__

#include <iostream>

#include "common.h"
#include "Composite.h"
#include "MCylinder.h"
#include "CCylinder.h"

///FORWARD DECLARATIONS
class Filament;
class Compartment;

///Cylinder class is a wrapper for a mechanical cylinder and chemical cylinder
/*!
 * Cylinder class is used to create a chemical and mechanical cylinder when needed.
 * It contains a constructor as well as getters for mcylinder and ccylinders.
 */

class Cylinder : public Composite {
    
private:
    Bead* _b1;  ///< Pointer to the first bead, associated with this cylinder ;
    Bead* _b2; ///< Pointer to the end bead in the cylinder.
               ///< Either empty - last cylinder, or pointer to the first Bead in a next cylinder.
    
    unique_ptr<MCylinder> _mCylinder; ///< ptr to mcylinder
    unique_ptr<CCylinder> _cCylinder; ///< ptr to ccylinder
    
    Filament* _pFilament; //< Pointer to filament where this cylinder belongs;
    int _positionFilament; ///< position on filament (1st, 2nd, ... etc)
    bool _last = false; ///< if the cylinder is last in the filament's cylinder list
    
    int _ID; ///Unique ID of cylinder, managed by CylinderDB
    
    Compartment* _compartment; ///< compartment this cylinder is currently in
    
    ///Function to find nearby cylinders in the grid, used in updatePosition
    vector<Cylinder*> findNearbyCylinders();
    
public:
    vector<double> coordinate; ///< coordinates of midpoint of cylinder, updated with updatePosition()
    
    ///Constructor and destructor
    Cylinder(Filament* f, Bead* b1, Bead* b2, int ID, bool extensionFront, bool extensionBack, bool creation);
    ~Cylinder();
    
    ///get mCylinder
    MCylinder* getMCylinder() {return _mCylinder.get();}
    
    ///get cCylinder
    CCylinder* getCCylinder() {return _cCylinder.get();}
    ///set cCylinder
    ///@note: since this is a unique ptr, will implicitly delete old CCylinder
    void setCCylinder(CCylinder* c) {_cCylinder = unique_ptr<CCylinder>(c);}
    
    ///get parent filament
    Filament* getFilament() {return _pFilament;}
    ///set parent filament
    void setFilament(Filament* pf) {_pFilament = pf;}
    
    ///Get beads
    Bead* getFirstBead() {return _b1;}
    Bead* getSecondBead() {return _b2;}
    
    Compartment* getCompartment() {return _compartment;}
    
    const int getID() {return _ID;}
    
    bool last(){ return _last;}
    void setLast(bool b){ _last = b;}
    
    void setPositionFilament(int positionFilament) {_positionFilament = positionFilament;}
    int getPositionFilament() {return _positionFilament;}
    
    ///Update the position of this cylinder
    ///@note - changes compartment of ccylinder if needed
    void updatePosition();
    
};



#endif /* defined(__Cyto__Cylinder__) */
