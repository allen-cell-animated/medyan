//
//  FilamentBending.h
//  Cyto
//
//  Created by Konstantin Popov on 8/27/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#ifndef Cyto_FilamentBending_h
#define Cyto_FilamentBending_h
#include "Bead.h"
#include "FilamentInteractions.h"
#include "Filament.h"

class FilamentInteractions;
class Filament;

template <class FBendingInteractionType>
class FilamentBending : public FilamentInteractions {
    
private:
    FBendingInteractionType _FFType;
    
public:
    
    double ComputeEnergy( Filament*, double d);
    void ComputeForces( Filament*);
    void ComputeForcesAux( Filament*);
};

#endif /* defined(__Cyto__FilamentBending__) */
