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
#include "ForceField.h"
#include <stdlib.h>

class BoundaryFF : public ForceField {
    
private:
    std::vector<std::unique_ptr<BoundaryInteraction>> _BoundaryInteractionVector;
    
public:
    BoundaryFF(std::string Interacion1, std::string Interacion2, std::string Interacion3 );
    
    // Public interfaecs to compute forces:
    virtual double ComputeEnergy(double d);
    virtual void ComputeForces();
    virtual void ComputeForcesAux();
};




#endif /* defined(__Cyto__BoundaryFF__) */
