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

class BoundaryElement;

class BoundaryInteractions {
private:
    std::string _name;
    
public:
    virtual double ComputeEnergy( BoundaryElement*,  double d) = 0;
    virtual void ComputeForces(BoundaryElement*) = 0;
    virtual void ComputeForcesAux(BoundaryElement*) = 0;
    
    // std::string getName() {return _name;}
    std::string getName() {return _name;}
    
};


#endif
