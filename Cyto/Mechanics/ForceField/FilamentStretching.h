//
//  FilamentStretching.h
//  Cyto
//
//  Created by Konstantin Popov on 8/19/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#ifndef Cyto_FilamentStretching_h
#define Cyto_FilamentStretching_h

#include "common.h"
#include "FilamentInteractions.h"


template <class FStretchingInteractionType>
class FilamentStretching : public FilamentInteractions {
    
private:
    FStretchingInteractionType _FFType;
    
public:
    virtual double computeEnergy( Filament*, double d);
    virtual void computeForces( Filament*);
    virtual void computeForcesAux( Filament*);
};


#endif /* defined(__Cyto__FilamentStretching__) */
