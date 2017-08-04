
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v3.1
//
//  Copyright (2015-2016)  Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the MEDYAN web page for more information:
//  http://www.medyan.org
//------------------------------------------------------------------

#ifndef MEDYAN_CGMethod_h
#define MEDYAN_CGMethod_h

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
    
    ///< Data vectors for calculation
    double *coord;  ///<bead coordinates (length 3*N)
    
    double *force; ///< bead forces (length 3*N)
    double *forceAux; ///< auxiliary force calculations (length 3*N)
    double *forceAuxPrev; ///<auxiliary force calculation previously (length 3*N)
    
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
    inline double allFDotF();
    inline double allFADotFA();
    inline double allFADotFAP();
    inline double allFDotFA();
    
    /// Get the max force in the system
    inline double maxF();
    
    /// Get bead with the max force in the system
    Bead* maxBead();
    
    /// Transfers data to lightweight arrays for min
    inline void startMinimization();
    /// Transfers updated coordinates and force to bead members
    inline void endMinimization();
    
    /// Move beads in search direction by d
    inline void moveBeads(double d);
    
    /// shift the gradient by d
    inline void shiftGradient(double d);
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
    
    /// Initialize data arrays
    inline void allocate(long numBeads) {
        
        coord = new double[3 * numBeads];
        force = new double[3 * numBeads];
        forceAux = new double[3 * numBeads];
        forceAuxPrev = new double[3 * numBeads];
    }
    
    ///Deallocation of CG arrays
    inline void deallocate() {
        
        delete coord;
        delete force;
        delete forceAux;
        delete forceAuxPrev;
    }
    
public:
    static long N; ///< Number of beads in the system, set before each minimization
    
    virtual ~CGMethod() {};
    
    /// Minimize the system
    virtual void minimize(ForceFieldManager &FFM, double GRADTOL,
                          double MAXDIST, double LAMBDAMAX, bool steplimit) = 0;
};


#endif
