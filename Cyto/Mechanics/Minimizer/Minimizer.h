//
//  Minimizer.h
//  Cyto
//
//  Created by Konstantin Popov on 9/8/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#ifndef __Cyto__Minimizer__
#define __Cyto__Minimizer__

#include <iostream>

#include "common.h"

///FORWARD DECLARATIONS
class ForceFieldManager;

class Minimizer {
public:
    void virtual equlibrate(ForceFieldManager &FFM) = 0;
};

#endif /* defined(__Cyto__Minimizer__) */
