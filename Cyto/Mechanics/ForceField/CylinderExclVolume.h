//
//  CylinderExclVolume.h
//  Cyto
//
//  Created by Konstantin Popov on 10/31/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#ifndef __Cyto__CylinderExclVolume__
#define __Cyto__CylinderExclVolume__

#include <stdio.h>
#include <iostream>

#include "common.h"
#include "CylinderVolumeInteractions.h"

class Cylinder;

template <class CVolumeInteractionType>
class CylinderExclVolume : public CylinderVolumeInteractions {
    
private:
    CVolumeInteractionType _FFType;
    
public:
    virtual double ComputeEnergy( Cylinder*, Cylinder*, double d);
    virtual void ComputeForces( Cylinder*, Cylinder*);
    virtual void ComputeForcesAux( Cylinder*, Cylinder* );
};



#endif /* defined(__Cyto__CylinderExclVolume__) */
