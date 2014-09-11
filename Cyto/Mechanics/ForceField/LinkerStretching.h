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
#include "Linker.h"
#include "LinkerInteractions.h"
#include "Bead.h"

class LinkerInteractions;
class Linker;

template <class LStretchingInteractionType>
class LinkerStretching : public LinkerInteractions
{
    
private:
    LStretchingInteractionType _FFType;
    
public:
    double ComputeEnergy( Linker*, double d);
    void ComputeForces( Linker*);
    void ComputeForcesAux( Linker*);
};

#endif /* defined(__Cyto__LinkerStratching__) */
