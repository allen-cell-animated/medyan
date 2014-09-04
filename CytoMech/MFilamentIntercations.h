//
//  MFilamentIntercations.h
//  Cyto
//
//  Created by Konstantin Popov on 8/15/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#ifndef __Cyto__MFilamentIntercations__
#define __Cyto__MFilamentIntercations__

#include <iostream>


class Filament;

class FilamentInteractions
{

    
private:
    
    std::string _name;
    
    
public:
    virtual double ComputeEnergy( Filament*,  double d) = 0;
    virtual double ComputeForces(Filament*) = 0;
    virtual double ComputeForcesAux(Filament*) = 0;
    
   // std::string getName() {return _name;}

    std::string getName() {return _name;}
    
};

#endif /* defined(__Cyto__MFilamentIntercations__) */
