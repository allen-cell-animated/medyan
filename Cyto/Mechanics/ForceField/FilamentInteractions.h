//
//  FilamentInteractions.h
//  Cyto
//
//  Created by Konstantin Popov on 8/15/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#ifndef __Cyto__FilamentInteractions__
#define __Cyto__FilamentInteractions__

#include <iostream>
#include "common.h"

class Filament;

class FilamentInteractions {
private:
    string _name;

public:
    virtual double computeEnergy( Filament*,  double d) = 0;
    virtual void computeForces(Filament*) = 0;
    virtual void computeForcesAux(Filament*) = 0;
    
   // string getName() {return _name;}
    const string& getName() {return _name;}
    
};

#endif /* defined(__Cyto__FilamentInteractions__) */
