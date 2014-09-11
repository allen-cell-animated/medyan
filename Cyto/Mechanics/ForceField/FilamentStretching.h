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

template <class FStretchingInteractionType>
class FilamentStretching : public FilamentInteractions {
    
private:
    FStretchingInteractionType _FFType;
    
public:
    virtual double ComputeEnergy( Filament*, double d);
    virtual void ComputeForces( Filament*);
    virtual void ComputeForcesAux( Filament*);
};


#endif /* defined(__Cyto__FilamentStretching__) */
