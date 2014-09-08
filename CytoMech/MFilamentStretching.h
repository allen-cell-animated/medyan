//
//  MFilamentStretching.h
//  Cyto
//
//  Created by Konstantin Popov on 8/19/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#ifndef Cyto_MFilamentStretching_h
#define Cyto_MFilamentStretching_h

#include "MFilament.h"
#include "MBead.h"
#include "MCylinder.h"


class FilamentInteractions;
class Filament;


template <class FStretchingInteractionType>
class FilamentStretching : private FilamentInteractions
{
    
private:
    FStretchingInteractionType _FFType;
    
    
public:
    
    double ComputeEnergy( Filament*, double d);
    
    void ComputeForces( Filament*);
    
    void ComputeForcesAux( Filament*);
    

};




#endif
