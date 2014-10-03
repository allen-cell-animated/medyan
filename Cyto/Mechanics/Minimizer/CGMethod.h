//
//  CGMethod.h
//  Cyto
//
//  Created by James Komianos on 9/10/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#ifndef __Cyto__CGMethod__
#define __Cyto__CGMethod__
#include <iostream>
#include <cmath>

#include "common.h"
#include "BeadDB.h"

class ForceFieldManager;

class CGMethod {
    
protected:
    double GradSquare();
    double GradSquare(int i);
    double GradDotProduct();
    void MoveBeads(double d);
    void ShiftGradient(double d);
    void PrintForces();
    double GoldenSection(ForceFieldManager &FFM);
    double GoldenSectionAlt(ForceFieldManager &FFM, double a, double b, double c, double tau);
    
    BeadDBKey getBeadDBKey() {return BeadDBKey();}
    
public:
    virtual void Minimize(ForceFieldManager &FFM) = 0;
};


#endif /* defined(__Cyto__CGMethod__) */
