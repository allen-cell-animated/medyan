//
//  MCGMethod.h
//  Cyto
//
//  Created by James Komianos on 9/10/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#ifndef __Cyto__MCGMethod__
#define __Cyto__MCGMethod__
#include <iostream>
#include "MBeadDB.h"
#include "ForceFieldManager.h"


class CGMethod {
    
protected:
    double GradSquare();
    double GradSquare(int i);
    double GradDotProduct();
    void MoveBeads(double d);
    void ShiftGradient(double d);
    void PrintForces();
    double GoldenSection(ForceFieldManager &FFM);
    
    BeadDBKey getBeadDBKey() {return BeadDBKey();}
    
public:
    virtual void Minimize(ForceFieldManager &FFM) = 0;
};


#endif /* defined(__Cyto__MCGMethod__) */
