//
//  MFilamentBending.h
//  Cyto
//
//  Created by Konstantin Popov on 8/27/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#ifndef Cyto_MFilamentBending_h
#define Cyto_MFilamentBending_h
#include "MBead.h"

class FilamentInteractions;
class Filament;


template <class BendingInteractionType>
class FilamentBending : private FilamentInteractions
{
    
private:
    BendingInteractionType _FFType;
    
    
public:
    
    double ComputeEnergy( Filament*, double d);
    
    void ComputeForces( Filament*);
    
    void ComputeForcesAux( Filament*);
    
    
};




#endif
