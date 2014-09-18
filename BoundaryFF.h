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

class BoundaryInteractions;

class BoundaryFF : public ForceField {
    
private:
    std::vector<std::unique_ptr<BoundaryInteractions>> _BoundaryInteractionVector;
    
public:
    BoundaryFF( std::string Interaction1, std::string Interaction2, std::string Interaction3 );
    
    // Public interfaecs to compute forces:
    virtual double ComputeEnergy(double d);
    virtual void ComputeForces();
    virtual void ComputeForcesAux();
};

#endif /* defined(__Cyto__BoundaryFF__) */
