//
//  CGPolakRibiereMethod.h
//  Cyto
//
//  Created by Konstantin Popov on 9/8/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#ifndef __Cyto__CGPolakRibiereMethod__
#define __Cyto__CGPolakRibiereMethod__

#include <iostream>

#include <iostream>
#include <cmath>
#include <numeric>
#include <vector>
#include <math.h>
#include <algorithm>

#include "common.h"

#include "CGMethod.h"

class PolakRibiere : public CGMethod
{
public:
    void minimize(ForceFieldManager &FFM);
};

#endif /* defined(__Cyto__CGPolakRibiereMethod__) */

