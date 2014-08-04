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
#include "Composite.h"

using namespace chem;

class Cylinder : public Composite {
    
private:
    MCylinder* _mCylinder; ///< ptr to mcylinder
    CCylinder* _cCylinder; ///< ptr to ccylinder
    
    Filament* _pFilament; //< Pointer to filament where this cylinder belongs;
    int _positionFilament; // position on filament (1st, 2nd, ... etc);
    
public:
    ///Constructor and destructor
    Cylinder(Filament* pf, Bead* firstBead, bool extension = false);
    ~Cylinder();
    
    ///get mCylinder
    MCylinder* getMCylinder() {return _mCylinder;}
    
    ///get cCylinder
    CCylinder* getCCylinder() {return _cCylinder;}
    
    bool IfLast();
    void SetLast(bool);
    
};



#endif /* defined(__Cyto__Cylinder__) */
