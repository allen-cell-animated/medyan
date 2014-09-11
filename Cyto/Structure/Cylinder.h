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
#include "MCylinder.h"
#include "CCylinder.h"
#include "ChemInitializer.h"
#include "Composite.h"

class Filament;

///Cylinder class is a wrapper for a mechanical cylinder and chemical cylinder
/*!
 * Cylinder class is used to create a chemical and mechanical cylinder when needed.
 * It contains a constructor as well as getters for mcylinder and ccylinders.
 */

class Cylinder : public Composite {
    
private:
    
    std::unique_ptr<MCylinder> _mCylinder; ///< ptr to mcylinder
    std::unique_ptr<CCylinder> _cCylinder; ///< ptr to ccylinder
    
    Filament* _pFilament; //< Pointer to filament where this cylinder belongs;
    int _positionFilament; // position on filament (1st, 2nd, ... etc);
    bool _ifLast;
    
public:
    ///Constructor and destructor
    Cylinder(Filament* pf, Bead* firstBead, Compartment* c, bool extensionFront, bool extensionBack);
    ~Cylinder();
    
    ///get mCylinder
    MCylinder* getMCylinder() {return _mCylinder.get();}
    
    ///get cCylinder
    CCylinder* getCCylinder() {return _cCylinder.get();}
    
    ///get parent filament
    Filament* getFilament() {return _pFilament;}
    
    bool IfLast();
    void SetLast(bool);
    
};



#endif /* defined(__Cyto__Cylinder__) */
