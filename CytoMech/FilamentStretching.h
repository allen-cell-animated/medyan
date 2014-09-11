//
//  FilamentStretching.h
//  Cyto
//
//  Created by Konstantin Popov on 8/19/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#ifndef Cyto_FilamentStretching_h
#define Cyto_FilamentStretching_h

#include "Filament.h"
#include "FilamentInteractions.h"
#include "Bead.h"
#include "Cylinder.h"

template <class FStretchingInteractionType>
class FilamentStretching : public FilamentInteractions {
    
private:
    FStretchingInteractionType _FFType;
    
public:
    double ComputeEnergy( Filament*, double d);
    void ComputeForces( Filament*);
    void ComputeForcesAux( Filament*);
};


#endif /* defined(__Cyto__FilamentStretching__) */
