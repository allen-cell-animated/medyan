//
//  MCGFletcherRievesMethod.h
//  Cyto
//
//  Created by Konstantin Popov on 9/8/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#ifndef __Cyto__MCGFletcherRievesMethod__
#define __Cyto__MCGFletcherRievesMethod__

#include <iostream>
#include <cmath>
#include <numeric>
#include <vector>
#include <math.h>
#include <algorithm>
#include "MCGMethod.h"

class FletcherRieves : public CGMethod
{
public:
   void Minimize(ForceFieldManager &FFM);
};

#endif /* defined(__Cyto__MCGFletcherRievesMethod__) */
