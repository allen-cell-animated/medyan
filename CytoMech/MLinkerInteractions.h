//
//  MLinkerInteractions.h
//  Cyto
//
//  Created by Konstantin Popov on 8/28/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#ifndef __Cyto__MLinkerInteractions__
#define __Cyto__MLinkerInteractions__

#include <iostream>



class Linker;

class LinkerInteractions
{
    
    
private:
    
    std::string _name;
    
    
public:
    virtual double ComputeEnergy( Linker*,  double d) = 0;
    virtual double ComputeForces(Linker*) = 0;
    virtual double ComputeForcesAux(Linker*) = 0;
    
    // std::string getName() {return _name;}
    
    std::string getName() {return _name;}
    
};














#endif /* defined(__Cyto__MLinkerInteractions__) */
