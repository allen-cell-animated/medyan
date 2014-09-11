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
#include "MotorGhost.h"
#include "MotorGhostInteractions.h"
#include "Bead.h"

class MotorGhostInteractions;

template <class MStretchingInteractionType>
class MotorGhostStretching : public MotorGhostInteractions {
    
private:
    MStretchingInteractionType _FFType;

public:
    double ComputeEnergy( MotorGhost*, double d);
    void ComputeForces( MotorGhost*);
    void ComputeForcesAux( MotorGhost*);
};



#endif /* defined(__Cyto__MotorGhostStretching__) */
