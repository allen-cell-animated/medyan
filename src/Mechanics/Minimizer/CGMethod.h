
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v3.2.1
//
//  Copyright (2015-2018)  Papoian Lab, University of Maryland
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
#ifdef CUDAACCL
#include "CUDAcommon.h"
#endif

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
    [[deprecated]] double *coord;  ///<bead coordinates (length 3*N)
    
    double *force = nullptr; ///< bead forces (length 3*N)
    double *forceAux = nullptr; ///< auxiliary force calculations (length 3*N)
    double *forceAuxPrev = nullptr; ///<auxiliary force calculation previously (length 3*N)
//    cylinder* cylindervec;

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

#ifdef CUDAACCL
    int *gpu_mutexlock;
    vector<int> blocksnthreads;
    vector<int> bntaddvector;
    vector<int> bnt;
    int *gpu_nint;
    double *gpu_g, *gpu_maxF;
    double *gSum;
    double *gSum2;
    double *gpu_fmax;
    double *g_currentenergy;
    double *gpu_params = NULL;
    double *gpu_FDotF;//curGrad
    double *gpu_FADotFA;//newGrad
    double *gpu_FADotFAP;//prevGrad
    double *gpu_FDotFA;
    double *gpu_initlambdalocal;
    bool *gpu_convergencecheck;
    bool *convergencecheck;
    double gpuFDotF(double *f1, double *f2);
    void CUDAresetlambda(cudaStream_t stream);
    void CUDAinitializeLambda(cudaStream_t stream1, bool *check_in, bool *check_out, bool
            *Polaksafestate, int *gpu_state);
    void CUDAfindLambda(cudaStream_t  stream1, cudaStream_t stream2, cudaEvent_t event, bool *checkin, bool
            *checkout, bool *gpu_safestate, int *gpu_state);
    void CUDAprepforbacktracking(cudaStream_t stream, bool *check_in, bool *check_out);
    void CUDAprepforsafebacktracking(cudaStream_t stream, bool *check_in, bool *check_out);
    void CUDAallFDotF(cudaStream_t stream);
    void CUDAallFADotFA(cudaStream_t stream);
    void CUDAallFADotFAP(cudaStream_t stream);
    void CUDAallFDotFA(cudaStream_t stream);
    void CUDAshiftGradient(cudaStream_t stream, bool *Mcheckin);
    void CUDAshiftGradientifSafe(cudaStream_t stream, bool *Mcheckin, bool *Scheckin);
//    void CUDAgetPolakvars(bool calc_safestate,cudaStream_t streamcalc, double* gpu_GRADTOL, bool *gminstatein,
//                                    bool *gminstateout, bool *gsafestateout, volatile bool *cminstate);
    void CUDAgetPolakvars(cudaStream_t streamcalc, double* gpu_GRADTOL, bool *gminstatein,
    bool *gminstateout, volatile bool *cminstate);
    void CUDAgetPolakvars2(cudaStream_t streamcalc, bool *gsafestateout);
    void CUDAinitializePolak(cudaStream_t stream, bool *minstatein, bool *minstateout, bool *safestatein, bool
    *safestateout);
    void CUDAmoveBeads(cudaStream_t stream, bool *gpu_checkin );
//    void getmaxFCUDA(double *gpu_forceAux, int *gpu_nint, double *gpu_fmax);
    //PING PONG for backtracking (both normal and safe)
//    struct backtrackingvars {
//        double currentEnergy;
//        double energyLambda;
//        double lambda;
//    };
    bool *g_stop1, *g_stop2, *g_s1, *g_s2, *g_ss;
//    backtrackingvars *bvar, *gpu_bvar1, *gpu_bvar2, *g_b1, *g_b2, *g_bs;
    cudaStream_t s1 = NULL, s2 = NULL, s3 = NULL, *sp1, *sp2, *sps, stream_bt = NULL;
    cudaStream_t stream_startmin = NULL;
    cudaEvent_t e1 = NULL, e2 = NULL, *ep1, *ep2, *eps;
    cudaEvent_t  e = NULL;
    int *gpu_state;
    // @PING PONG ENDS

#endif
    volatile bool *cconvergencecheck;
    bool *h_stop, sconvergencecheck;
    /// For use in minimization


    double allFDotF();
    double allFADotFA();
    double allFADotFAP();
    double allFDotFA();
    
    /// Get the max force in the system
    double maxF();
    
    /// Get bead with the max force in the system
    Bead* maxBead();
    
    /// Transfers data to lightweight arrays for min
    void startMinimization();
    /// Transfers updated coordinates and force to bead members
    void endMinimization();
    
    /// Move beads in search direction by d
    void moveBeads(double d);
    
    /// shift the gradient by d
    void shiftGradient(double d);
    //@}

#ifdef CUDAACCL
    //@{
    double backtrackingLineSearchCUDA(ForceFieldManager& FFM, double MAXDIST,
                                  double LAMBDAMAX, bool *gpu_safestate);
    //@}
#endif // CUDAACCL

    //@{
    /// Linear search methods
    /// A simple backtracking search method that computes an optimal
    /// energy change and compares the backtracked energy to it
    double backtrackingLineSearch(ForceFieldManager& FFM, double MAXDIST,
                                                          double LAMBDAMAX, bool *gpu_safestate);
    
    /// The safemode backtracking search, returns the first energy decrease
    ///@note - The most robust linesearch method, but very slow

    double safeBacktrackingLineSearch(ForceFieldManager& FFM, double MAXDIST,
                                                              double LAMBDAMAX, bool *gpu_safestate);

    //@}
    
    /// Print forces on all beads
    void printForces();
    
public:
    static long N; ///< Number of beads in the system, set before each minimization
    static long Ncyl;
    
    virtual ~CGMethod() {};
    
    /// Minimize the system
    virtual void minimize(ForceFieldManager &FFM, double GRADTOL,
                          double MAXDIST, double LAMBDAMAX, bool steplimit) = 0;
};


#endif
