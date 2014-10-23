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

class Filament;
class Compartment;

///Cylinder class is a wrapper for a mechanical cylinder and chemical cylinder
/*!
 * Cylinder class is used to create a chemical and mechanical cylinder when needed.
 * It contains a constructor as well as getters for mcylinder and ccylinders.
 */

class Cylinder : public Composite {
    
private:
    Bead* _pFirst;  ///< Pointer to the first bead, associated with this cylinder ;
    Bead* _pSecond; ///< Pointer to the end bead in the cylinder. Either empty - last cylinder, or pointer to the first Bead in a next cylinder.
    
    std::unique_ptr<MCylinder> _mCylinder; ///< ptr to mcylinder
    std::unique_ptr<CCylinder> _cCylinder; ///< ptr to ccylinder
    
    Filament* _pFilament; //< Pointer to filament where this cylinder belongs;
    int _positionFilament; ///< position on filament (1st, 2nd, ... etc)
    bool _ifLast = false; ///< if the cylinder is last in the filament's cylinder list
    
    Compartment* _compartment; ///< compartment this cylinder is currently in
    
public:
    ///Constructor and destructor
    Cylinder(Filament* pf, Bead* firstBead, Bead* secondBead, bool extensionFront, bool extensionBack, bool creation);
    ~Cylinder();
    
    ///get mCylinder
    MCylinder* getMCylinder() {return _mCylinder.get();}
    
    ///get cCylinder
    CCylinder* getCCylinder() {return _cCylinder.get();}
    ///set cCylinder
    ///@note: since this is a unique ptr, will implicitly delete old CCylinder
    void setCCylinder(CCylinder* c) {_cCylinder = std::unique_ptr<CCylinder>(c);}
    
    ///get parent filament
    Filament* getFilament() {return _pFilament;}
    ///set parent filament
    void setFilament(Filament* pf) {_pFilament = pf;}
    
    ///Get beads
    Bead* GetFirstBead() {return _pFirst;}
    Bead* GetSecondBead() {return _pSecond;}
    
    bool IfLast();
    void SetLast(bool);
    
    ///Update the position of this cylinder
    ///@note - changes compartment of ccylinder if needed
    void updatePosition();
    
};



#endif /* defined(__Cyto__Cylinder__) */
