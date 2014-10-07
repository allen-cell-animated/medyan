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

#include "BeadDB.h"
#include "common.h"

class ForceFieldManager;

class CGMethod {
    
protected:
    ///helpers for searching and bracketing
    void swap(double &a, double &b);
    void shift2(double &a, double &b, double c);
    void shift3(double &a, double &b, double &c, double d);
    double sign(double a, double b);
    
    ///bracketing function (from Numerical Recipes in C++, second edition)
    void makeBracket(ForceFieldManager &FFM, double &ax, double &bx, double &cx, double &fa, double &fb, double &fc);
    
    ///Gradient calculations
    double GradSquare();
    double GradSquare(int i);
    double GradDotProduct();
    void MoveBeads(double d);
    void ShiftGradient(double d);
    
    ///various linear search methods
    double GoldenSection1(ForceFieldManager &FFM);
    double GoldenSection2(ForceFieldManager &FFM, double ax, double bx, double cx, double tau);
    double GoldenSection3(ForceFieldManager &FFM, double ax, double bx, double cx, double tol);
    
    double BinarySearch(ForceFieldManager& FFM, double a, double b );
    
    BeadDBKey getBeadDBKey() {return BeadDBKey();}
    
    void PrintForces();
    
public:
    virtual void Minimize(ForceFieldManager &FFM) = 0;
};


#endif /* defined(__Cyto__CGMethod__) */
