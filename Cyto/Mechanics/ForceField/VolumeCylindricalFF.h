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

class CylinderVolumeInteractions;

class VolumeCylindricalFF : public ForceField {
    
private:
    std::vector <std::unique_ptr<CylinderVolumeInteractions>> _cylinderVolInteractionVector;
    
public:
    VolumeCylindricalFF(std::string& Interaction );
    
    // Public interfaecs to compute forces:
    virtual double ComputeEnergy(double d);
    virtual void ComputeForces();
    virtual void ComputeForcesAux();
};




#endif /* defined(__Cyto__VolumeCylindricalFF__) */
