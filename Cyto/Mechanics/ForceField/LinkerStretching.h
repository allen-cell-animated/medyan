//
//  LinkerStretching.h
//  Cyto
//
//  Created by Konstantin Popov on 8/28/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#ifndef __Cyto__LinkerStratching__
#define __Cyto__LinkerStratching__

#include <iostream>

#include "common.h"
#include "LinkerInteractions.h"

///FORWARD DECLARATIONS
class Linker;

template <class LStretchingInteractionType>
class LinkerStretching : public LinkerInteractions {
    
private:
    LStretchingInteractionType _FFType;
    
public:
    virtual double computeEnergy(Linker*, double d);
    virtual void computeForces(Linker*);
    virtual void computeForcesAux(Linker*);
};

#endif /* defined(__Cyto__LinkerStretching__) */
