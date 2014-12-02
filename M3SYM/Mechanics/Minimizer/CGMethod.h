
//------------------------------------------------------------------
//  **M3SYM** - Simulation Package for the Mechanochemical
//              Dynamics of Active Networks, 3rd Generation
//
//  Copyright (2014) Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the Papoian lab page for installation and documentation:
//  http://papoian.chem.umd.edu/
//------------------------------------------------------------------

#ifndef M3SYM_CGMethod_h
#define M3SYM_CGMethod_h

#include <cmath>

#include "common.h"

//FORWARD DECLARATIONS
class ForceFieldManager;

/// For performing a conjugate gradient minimization method

/*!
 *  CGMethod is an abstract class that contains methods for conjugate gradient minimization. It has functions
 *  for various line search techniques, including golden section, backtracking, and quadratic line search.
 *  This base class also contains parameters for tolerance criterion and other helpful parameters.
 */
class CGMethod {

protected:
    
    //@{
    /// Lambda parameter for use in linear search methods
    const double LAMBDAMIN = 0.001; ///< Minimum lambda that can be returned, used only in golden section for now
    const double LAMBDAMAX = 1; ///< Max lambda that can be returned, used in all methods
    const double MAXDIST = 1; ///< Max distance beads can be moved, used only in backtracking line search
    //@}
    
    //@{
    /// Parameter used in backtracking line search
    const double LAMBDAREDUCE = 0.1; ///< Lambda reduction parameter for backtracking
    const double BACKTRACKSLOPE = 0.1; ///< Backtrack slope parameter
    //@}
    
    //@{
    /// Parameter used in quadratic line search
    const double QUADRATICTOL = 0.1;
    const double EPSQUAD = 1e-28;
    //@}
    
    //@{
    /// Parameter used in golden section
    const double PHI = (1 + sqrt(5)) / 2;
    const double R = 0.61803399;
    const double C = 1 - R;
    //@}
    
    const double LSENERGYTOL = 1e-6; ///< Line search energy tolerance for all linesearch methods
    
    ///helpers for searching and bracketing
    void swap(double &a, double &b);
    void shift2(double &a, double &b, double c);
    void shift3(double &a, double &b, double &c, double d);
    double sign(double a, double b);
    
protected:
    
    ///Energy counter, for use in linear search methods
    int _energyChangeCounter = 0; ///< Number of iterations where energy has not changed by an amount more than LSENERGYTOL
    const int ENERGYCHANGEITER = 20; ///< Max number of iterations allowed where """
    
    const double GRADTOL = 1e-10; ///< Gradient minimization tolerance
    
    /// Gracketing function (from Numerical Recipes in C++, second edition)
    void makeBracket(ForceFieldManager &FFM, double &ax, double &bx, double &cx, double &fa, double &fb, double &fc);
    
    //@{
    ///For use in minimization
    double gradSquare();
    double gradAuxSquare();
    double gradDotProduct();
    
    void moveBeads(double d);
    void moveBeadsAux(double d);
    void shiftGradient(double d);
    //@}
    
    //@{
    /// Linear search method
    double goldenSection1(ForceFieldManager &FFM);
    double goldenSection2(ForceFieldManager &FFM);
    double binarySearch(ForceFieldManager& FFM);
    
    double backtrackingLineSearch(ForceFieldManager& FFM);
    double quadraticLineSearch(ForceFieldManager& FFM);
    //@}
    
    void printForces();
    
public:
    /// Minimize the system
    virtual void minimize(ForceFieldManager &FFM) = 0;
};


#endif
