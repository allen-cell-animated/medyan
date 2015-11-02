
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
class Bead;

/// For performing a conjugate gradient minimization method

/*!
 *  CGMethod is an abstract class that contains methods for conjugate gradient 
 *  minimization. It has functions for various line search techniques, including golden 
 *  section, backtracking, and quadratic line search. This base class also contains 
 *  parameters for tolerance criterion and other helpful parameters.
 */
class CGMethod {

protected:
    
    /// Safe mode which chooses the safe backtracking search if the
    /// minimizer got itself into trouble.
    bool _safeMode = false;
    
    /// Minimum number of minimization steps, in the case of a
    /// small number of beads in the system
    const int _MINNUMSTEPS = 1E4;
    
    //@{
    /// Parameter used in backtracking line search
    const double LAMBDAREDUCE = 0.5;     ///< Lambda reduction parameter for backtracking
    const double LAMBDATOL = 1e-8;       ///< Lambda tolerance parameter
    
    const double SAFELAMBDAREDUCE = 0.9;  ///< Lambda reduction parameter for conservative backtracking
    
    const double BACKTRACKSLOPE = 0.4;   ///< Backtracking slope
    //@}
    
    //@{
    /// For use in minimization
    double allFDotF();
    double allFADotFA();
    double allFADotFAP();
    double allFDotFA();
    
    /// Get the max force in the system
    double maxF();
    
    /// Get bead with the max force in the system
    Bead* maxBead();
    
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
    /// A simple backtracking search method that computes an optimal
    /// energy change and compares the backtracked energy to it
    double backtrackingLineSearch(ForceFieldManager& FFM, double MAXDIST,
                                                          double LAMBDAMAX);
    
    /// The safemode backtracking search, returns the first energy decrease
    ///@note - The most robust linesearch method, but very slow
    double safeBacktrackingLineSearch(ForceFieldManager& FFM, double MAXDIST,
                                                              double LAMBDAMAX);
    //@}
    
    /// Print forces on all beads
    void printForces();
    
public:
    
    virtual ~CGMethod() {};
    
    /// Minimize the system
    virtual void minimize(ForceFieldManager &FFM, double GRADTOL,
                          double MAXDIST, double LAMBDAMAX, bool steplimit) = 0;
};


#endif
