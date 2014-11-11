//
//  VolumeCylindricalFF.h
//  Cyto
//
//  Created by Konstantin Popov on 10/29/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#ifndef __Cyto__VolumeCylindricalFF__
#define __Cyto__VolumeCylindricalFF__

#include <stdio.h>
#include <iostream>
#include <vector>

#include "common.h"
#include "ForceField.h"

///FORWARD DECLARATIONS
class CylinderVolumeInteractions;

class VolumeCylindricalFF : public ForceField {
    
private:
    vector <unique_ptr<CylinderVolumeInteractions>> _cylinderVolInteractionVector;
    
public:
    VolumeCylindricalFF(string& interaction );
    
    // Public interfaces to compute forces:
    virtual double computeEnergy(double d);
    virtual void computeForces();
    virtual void computeForcesAux();
};




#endif /* defined(__Cyto__VolumeCylindricalFF__) */
