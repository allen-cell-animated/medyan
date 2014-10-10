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


///CGMethod class is for performing a conjugate gradient method
/*!
 *  CGMethod is an abstract class that contains methods for conjugate gradient minimization. It has functions
 *  for various line search techniques, including golden section, backtracking, and quadratic line search.
 *  This base class also contains parameters for tolerance criterion and other helpful parameters.
 */

class CGMethod {

private:
    
    ///Lambda parameters for use in linear search methods
    const double LAMBDAMIN = 0.001; ///< Minimum lambda that can be returned, used only in golden section for now
    const double LAMBDAMAX = 10; ///< Max lambda that can be returned, used in all methods
    const double MAXDIST = 10; ///< Max distance beads can be moved, used only in backtracking line search
    
    ///Parameters used in backtracking line search
    const double LAMBDAREDUCE = 0.1; ///< Lambda reduction parameter for backtracking
    const double BACKTRACKSLOPE = 0.1; ///< Backtrack slope parameter
    
    ///Parameters used in quadratic line search
    const double QUADRATICTOL = 0.1;
    const double EPSQUAD = 1e-28;
    
    ///Parameters used in golden section
    const double PHI = (1 + sqrt(5)) / 2;
    const double R = 0.61803399;
    const double C = 1 - R;
    
    const double LSENERGYTOL = 1e-4; ///<Line search energy tolerance for all linesearch methods
    
    ///helpers for searching and bracketing
    void swap(double &a, double &b);
    void shift2(double &a, double &b, double c);
    void shift3(double &a, double &b, double &c, double d);
    double sign(double a, double b);
    
protected:
    
    ///Energy counter, for use in linear search methods
    int _energyChangeCounter = 0; ///<number of iterations where energy has not changed by an amount more than LSENERGYTOL
    const int ENERGYCHANGEITER = 20; ///<max number of iterations allowed where """
    
    const double GRADTOL = 1e-4; ///< gradient minimization tolerance
    
    ///bracketing function (from Numerical Recipes in C++, second edition)
    void makeBracket(ForceFieldManager &FFM, double &ax, double &bx, double &cx, double &fa, double &fb, double &fc);
    
    ///Gradient calculations
    double GradSquare();
    double GradAuxSquare();
    double GradDotProduct();
    void MoveBeads(double d);
    void MoveBeadsAux(double d);
    void ShiftGradient(double d);
    
    ///various linear search methods
    double GoldenSection1(ForceFieldManager &FFM);
    double GoldenSection2(ForceFieldManager &FFM);
    double BinarySearch(ForceFieldManager& FFM);
    
    double BacktrackingLineSearch(ForceFieldManager& FFM);
    double QuadraticLineSearch(ForceFieldManager& FFM);
    
    BeadDBKey getBeadDBKey() {return BeadDBKey();}
    
    void PrintForces();
    
public:
    virtual void Minimize(ForceFieldManager &FFM) = 0;
};


#endif /* defined(__Cyto__CGMethod__) */
