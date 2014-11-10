//
//  ForceField.h
//  Cyto
//
//  Created by James Komianos on 8/5/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#ifndef __Cyto__ForceField__
#define __Cyto__ForceField__

#include <iostream>

#include "common.h"

///ForceField is an abstract class to represent various force field calculations
/*!
 *  ForceField is used for force calculations between filaments, beads, linkers, etc.
 *  Specific implementations of the ForceField class will have different potentials.
 */
class ForceField {

private:
    string _name;
    
public:
    const string& getName() {return _name;}
    virtual double computeEnergy(double d) = 0;
    virtual void computeForces() = 0;
    virtual void computeForcesAux() = 0;
};






#endif /* defined(__Cyto__ForceField__) */
