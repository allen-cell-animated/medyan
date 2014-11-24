//
//  MotorGhostStretching.h
//  Cyto
//
//  Created by Konstantin Popov on 9/3/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#ifndef __Cyto__MotorGhostStretching__
#define __Cyto__MotorGhostStretching__

#include <iostream>

#include "common.h"

#include "MotorGhostInteractions.h"

///FORWARD DECLARATIONS
class MotorGhost;

template <class MStretchingInteractionType>
class MotorGhostStretching : public MotorGhostInteractions {
    
private:
    MStretchingInteractionType _FFType;

public:
    virtual double computeEnergy(MotorGhost*, double d);
    virtual void computeForces(MotorGhost*);
    virtual void computeForcesAux(MotorGhost*);
};



#endif /* defined(__Cyto__MotorGhostStretching__) */
