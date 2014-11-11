//
//  BoundaryInteractions.h
//  Cyto
//
//  Created by Konstantin Popov on 9/12/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#ifndef Cyto_BoundaryInteractions_h
#define Cyto_BoundaryInteractions_h

#include "common.h"
#include <iostream>

///FORWARD DECLARATIONS
class BoundaryElement;

class BoundaryInteractions {
private:
    string _name;
    
public:
    virtual double computeEnergy(BoundaryElement*,  double d) = 0;
    virtual void computeForces(BoundaryElement*) = 0;
    virtual void computeForcesAux(BoundaryElement*) = 0;
    
    // string getName() {return _name;}
    string getName() {return _name;}
    
};


#endif
