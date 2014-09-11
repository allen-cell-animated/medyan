//
//  CGFletcherRievesMethod.h
//  Cyto
//
//  Created by Konstantin Popov on 9/8/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#ifndef __Cyto__CGFletcherRievesMethod__
#define __Cyto__CGFletcherRievesMethod__

#include <iostream>
#include <cmath>
#include <numeric>
#include <vector>
#include <math.h>
#include <algorithm>
#include "CGMethod.h"

class FletcherRieves : public CGMethod
{
public:
   void Minimize(ForceFieldManager &FFM);
};

#endif /* defined(__Cyto__CGFletcherRievesMethod__) */
