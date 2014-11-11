//
//  BoundaryFF.h
//  Cyto
//
//  Created by Konstantin Popov on 9/11/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#ifndef __Cyto__BoundaryFF__
#define __Cyto__BoundaryFF__

#include <iostream>
#include <vector>
#include <stdlib.h>

#include "common.h"
#include "ForceField.h"

///FORWARD DECLARATIONS
class BoundaryInteractions;

class BoundaryFF : public ForceField {
    
private:
    vector<unique_ptr<BoundaryInteractions>> _BoundaryInteractionVector;
    
public:
    BoundaryFF( string interaction1, string interaction2, string interaction3 );
    
    // Public interfaces to compute forces:
    virtual double computeEnergy(double d);
    virtual void computeForces();
    virtual void computeForcesAux();
};

#endif /* defined(__Cyto__BoundaryFF__) */
