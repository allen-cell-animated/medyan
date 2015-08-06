
//------------------------------------------------------------------
//  **M3SYM** - Simulation Package for the Mechanochemical
//              Dynamics of Active Networks, 3rd Generation
//
//  Copyright (2014) Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the Papoian lab page for more information:
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
 *  CGMethod is an abstract class that contains methods for conjugate gradient 
 *  minimization. It has functions for various line search techniques, including golden 
 *  section, backtracking, and quadratic line search. This base class also contains 
 *  parameters for tolerance criterion and other helpful parameters.
 */
class CGMethod {

protected:
    
    //@{
    /// Parameter used in backtracking line search
    const double LAMBDAREDUCE = 0.5;   ///< Lambda reduction parameter for backtracking
    const double LAMBDATOL = 1e-12;    ///< Lambda tolerance parameter
    //@}
    
    //@{
    /// For use in minimization
    double allFDotF();
    double allFADotFA();
    double allFADotFAP();
    double allFDotFA();
    
    /// Get the max force in the system
    double maxF();
    
    /// Sets coordinates before minimization
    void startMinimization();
    
    /// Move beads in search direction by d
    void moveBeads(double d);
    /// Reset to previous position
    void resetBeads();
    /// Update the previous position
    void setBeads();
    
    /// shift the gradient by d
    void shiftGradient(double d);
    //@}
    
    //@{
    /// Linear search methods
    ///@note - The most robust linesearch method, but slow at times.
    double backtrackingLineSearch(ForceFieldManager& FFM, double MAXDIST,
                                                          double LAMBDAMAX);
    //@}
    
    /// Print forces on all beads
    void printForces();
    
public:
    
    virtual ~CGMethod() {};
    
    /// Minimize the system
    virtual void minimize(ForceFieldManager &FFM, double GRADTOL,
                                                  double MAXDIST,
                                                  double LAMBDAMAX) = 0;
};


#endif
