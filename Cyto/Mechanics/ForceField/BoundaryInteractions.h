//
//  BoundaryInteractions.h
//  Cyto
//
//  Created by Konstantin Popov on 9/12/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#ifndef Cyto_BoundaryInteractions_h
#define Cyto_BoundaryInteractions_h

#include <iostream>

#include "common.h"

#include "NeighborListContainer.h"
#include "SystemParameters.h"

//FORWARD DECLARATIONS
class BoundaryElement;
class Bead;

class BoundaryInteractions : public BoundaryElementNLContainer {
private:
    string _name;
    
public:
    BoundaryInteractions() : BoundaryElementNLContainer(SystemParameters::Boundaries().boundaryCutoff) {}
    
    virtual double computeEnergy(BoundaryElement*, Bead*, double d) = 0;
    virtual void computeForces(BoundaryElement*, Bead*) = 0;
    virtual void computeForcesAux(BoundaryElement*, Bead*) = 0;
    
    // string getName() {return _name;}
    string getName() {return _name;}
    
};


#endif
