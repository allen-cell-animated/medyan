//
//  MCGPolakRibiereMethod.h
//  Cyto
//
//  Created by Konstantin Popov on 9/8/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#ifndef __Cyto__MCGPolakRibiereMethod__
#define __Cyto__MCGPolakRibiereMethod__

#include <iostream>

#endif /* defined(__Cyto__MCGPolakRibiereMethod__) */

#include <iostream>
#include <cmath>
#include <numeric>
#include <vector>
#include <math.h>
#include <algorithm>
#include "Mcommon.h"
#include "MFilament.h"
#include "SubSystem.h"

using namespace std;

class MController;

class PolakRibiere
{
    
private:
    double GradSquare();
    double GradSquare(int i);
    double GradDotProduct();
    void MoveBeads(double d);
    void ShiftGradient(double d);
    void PrintForces();
    double GoldenSection(MController* mc);
    
public:
    
    void Minimize(MController* mc);
};
