
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
#include "CUDAcommon.h"

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
    floatingpoint *coord;  ///<bead coordinates (length 3*N)

    floatingpoint *force; ///< bead forces (length 3*N)
    floatingpoint *forceAux; ///< auxiliary force calculations (length 3*N)
    floatingpoint *forceAuxPrev; ///<auxiliary force calculation previously (length 3*N)
//    cylinder* cylindervec;

    /// Safe mode which chooses the safe backtracking search if the
    /// minimizer got itself into trouble.
    bool _safeMode = false;
    
    /// Minimum number of minimization steps, in the case of a
    /// small number of beads in the system
    const int _MINNUMSTEPS = 1E4;
    
    //@{
    /// Parameter used in backtracking line search
    const floatingpoint LAMBDAREDUCE = 0.5;     ///< Lambda reduction parameter for backtracking
    floatingpoint LAMBDATOL = 1e-4;       ///< Lambda tolerance parameter
    
    const floatingpoint SAFELAMBDAREDUCE = 0.9;  ///< Lambda reduction parameter for conservative backtracking
    
    const floatingpoint BACKTRACKSLOPE = 0.4;   ///< Backtracking slope
    //@}
    
    //@{

#ifdef CUDAACCL
    int *gpu_mutexlock;
    vector<int> blocksnthreads;
    vector<int> bntaddvector;
    vector<int> bnt;
    int *gpu_nint;
    floatingpoint *gpu_g, *gpu_maxF;
    floatingpoint *gSum;
    floatingpoint *gSum2;
    floatingpoint *gpu_fmax;
    floatingpoint *g_currentenergy;
    floatingpoint *gpu_params = NULL;
    floatingpoint *gpu_FDotF;//curGrad
    floatingpoint *gpu_FADotFA;//newGrad
    floatingpoint *gpu_FADotFAP;//prevGrad
    floatingpoint *gpu_FDotFA;
    floatingpoint *gpu_initlambdalocal;
    bool *gpu_convergencecheck;
    bool *convergencecheck;
    floatingpoint gpuFDotF(floatingpoint *f1, floatingpoint *f2);
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
//    void CUDAgetPolakvars(bool calc_safestate,cudaStream_t streamcalc, floatingpoint* gpu_GRADTOL, bool *gminstatein,
//                                    bool *gminstateout, bool *gsafestateout, volatile bool *cminstate);
    void CUDAgetPolakvars(cudaStream_t streamcalc, floatingpoint* gpu_GRADTOL, bool *gminstatein,
    bool *gminstateout, volatile bool *cminstate);
    void CUDAgetPolakvars2(cudaStream_t streamcalc, bool *gsafestateout);
    void CUDAinitializePolak(cudaStream_t stream, bool *minstatein, bool *minstateout, bool *safestatein, bool
    *safestateout);
    void CUDAmoveBeads(cudaStream_t stream, bool *gpu_checkin );
//    void getmaxFCUDA(floatingpoint *gpu_forceAux, int *gpu_nint, floatingpoint *gpu_fmax);
    //PING PONG for backtracking (both normal and safe)
//    struct backtrackingvars {
//        floatingpoint currentEnergy;
//        floatingpoint energyLambda;
//        floatingpoint lambda;
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


    floatingpoint allFDotF();
    floatingpoint allFADotFA();
    floatingpoint allFADotFAP();
    floatingpoint allFDotFA();
    
    /// Get the max force in the system
    floatingpoint maxF();
    
    /// Get bead with the max force in the system
    Bead* maxBead();
    
    /// Transfers data to lightweight arrays for min
    void startMinimization();
    /// Transfers updated coordinates and force to bead members
    void endMinimization();
    
    /// Move beads in search direction by d
    void moveBeads(floatingpoint d);

    /// shift the gradient by d
    void shiftGradient(floatingpoint d);
    //@}

#ifdef CUDAACCL
    //@{
    floatingpoint backtrackingLineSearchCUDA(ForceFieldManager& FFM, floatingpoint MAXDIST,
                                  floatingpoint LAMBDAMAX, bool *gpu_safestate);
    //@}
#endif // CUDAACCL

    //@{
    /// Linear search methods
    /// A simple backtracking search method that computes an optimal
    /// energy change and compares the backtracked energy to it
    floatingpoint backtrackingLineSearch(ForceFieldManager& FFM, floatingpoint MAXDIST,
                                                          floatingpoint LAMBDAMAX, bool *gpu_safestate);
    
    /// The safemode backtracking search, returns the first energy decrease
    ///@note - The most robust linesearch method, but very slow

    floatingpoint safeBacktrackingLineSearch(ForceFieldManager& FFM,
    		floatingpoint MAXDIST, floatingpoint LAMBDAMAX, bool *gpu_safestate);

    floatingpoint setLAMBDATOL(int maxF_order){

        int orderdimension = 3; ///1000s of nm
        int LAMBDATOLorder = -(6-orderdimension) - maxF_order;
        LAMBDATOL = 1;
        if(LAMBDATOLorder > 0){
            for(int i =0; i < LAMBDATOLorder; i ++)
                LAMBDATOL *= 10;
        }
        else{
            for(int i =0; i > LAMBDATOLorder; i --)
                LAMBDATOL *= 0.1;
        }

        LAMBDATOL = max<floatingpoint>(1e-8, LAMBDATOL);
        LAMBDATOL = min<floatingpoint>(1e-1, LAMBDATOL);

        cout<<"maxF order "<<maxF_order<<" lambdatol "<<LAMBDATOL<<endl;
    }

    //@}
    
    /// Print forces on all beads
    void printForces();
    
    /// Initialize data arrays
    inline void allocate(long numBeadsx3, long Ncyl) {

//        coord = new floatingpoint[numBeadsx3];
        force = new floatingpoint[numBeadsx3];
        forceAux = new floatingpoint[numBeadsx3];
        forceAuxPrev = new floatingpoint[numBeadsx3];

        for(int i =0; i < numBeadsx3; i++){
        	force[i] = 0.0;
        	forceAux[i]=0.0;
        	forceAuxPrev[i]=0.0;
        }
    }
    
    ///Deallocation of CG arrays
    inline void deallocate() {
//        coord = CUDAcommon::serlvars.coord;
//        delete [] coord;
        delete [] force;
        delete [] forceAux;
        delete [] forceAuxPrev;
    }
public:
    static long N; ///< Number of beads in the system, set before each minimization
    static long Ncyl;
    
    virtual ~CGMethod() {};
    
    /// Minimize the system
    virtual void minimize(ForceFieldManager &FFM, floatingpoint GRADTOL,
                          floatingpoint MAXDIST, floatingpoint LAMBDAMAX, bool steplimit) = 0;
};


#endif
