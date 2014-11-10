//
//  LinkerFF.h
//  Cyto
//
//  Created by Konstantin Popov on 9/5/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#ifndef __Cyto__LinkerFF__
#define __Cyto__LinkerFF__

#include <iostream>
#include <vector>

#include "common.h"
#include "ForceField.h"

class LinkerInteractions;

class LinkerFF : public ForceField {
    
private:
    vector<unique_ptr<LinkerInteractions>> _linkerInteractionVector;
    
public:
    LinkerFF(string& stretching, string& bending, string& twisting );
    
    // Public interfaces to compute forces:
    virtual double computeEnergy(double d);
    virtual void computeForces();
    virtual void computeForcesAux();
    
};


#endif /* defined(__Cyto__LinkerFF__) */
