//
//  LinkerInteractions.h
//  Cyto
//
//  Created by Konstantin Popov on 8/28/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#ifndef __Cyto__LinkerInteractions__
#define __Cyto__LinkerInteractions__

#include <iostream>

#include "common.h"

///FORWARD DECLARATIONS
class Linker;

class LinkerInteractions {
private:
    string _name;
    
public:
    virtual double computeEnergy(Linker*,  double d) = 0;
    virtual void computeForces(Linker*) = 0;
    virtual void computeForcesAux(Linker*) = 0;
    
    // string getName() {return _name;}
    
    const string& getName() {return _name;}
    
};

#endif /* defined(__Cyto__LinkerInteractions__) */
