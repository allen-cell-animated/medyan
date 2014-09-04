//
//  MMotorGhostStretching.h
//  Cyto
//
//  Created by Konstantin Popov on 9/3/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#ifndef __Cyto__MMotorGhostStretching__
#define __Cyto__MMotorGhostStretching__

#include <iostream>
#include "MMotorGhost.h"
#include "MBead.h"


class LinkerInteractions;
class Linker;


template <class StretchingInteractionType>
class MotorGhostStretching : private MotorGhostStretching
{
    
private:
    StretchingInteractionType _FFType;
    
    
public:
    
    double ComputeEnergy( MotorGhost*, double d);
    
    void ComputeForces( MotorGhost*);
    
    void ComputeForcesAux( MotorGhost*);
    
    
};



#endif /* defined(__Cyto__MMotorGhostStretching__) */
