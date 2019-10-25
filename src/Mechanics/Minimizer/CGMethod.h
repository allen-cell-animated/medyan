
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v4.0
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
#include "CUDAcommon.h"
#include "Mechanics/Minimizer/MinimizationTypes.hpp"

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
    chrono::high_resolution_clock::time_point tbegin, tend;

    ///< Data vectors for calculation
    [[deprecated]] floatingpoint *coord;  ///<bead coordinates (length 3*N)
    [[deprecated]] floatingpoint *coordlineSearch; ///coords used during line search

    [[deprecated]] floatingpoint *force=NULL; ///< bead forces (length 3*N)
    [[deprecated]] floatingpoint *forceAux=NULL; ///< auxiliary force calculations (length 3*N)
    [[deprecated]] floatingpoint *forceAuxPrev=NULL; ///<auxiliary force calculation previously (length
    // 3*N)
//    cylinder* cylindervec;

//Gradients
	floatingpoint FADotFA = 0.0;
	floatingpoint FADotFAP = 0.0;
    /// Safe mode which chooses the safe backtracking search if the
    /// minimizer got itself into trouble.
    bool _safeMode = false;

    /// Minimum number of minimization steps, in the case of a
    /// small number of beads in the system
    const int _MINNUMSTEPS = 1E4;

    //@{
    /// Parameter used in backtracking line search
    const floatingpoint LAMBDAREDUCE = 0.5;     ///< Lambda reduction parameter for backtracking
    floatingpoint LAMBDATOL = 1e-8;       ///< Lambda tolerance parameter

    const floatingpoint SAFELAMBDAREDUCE = 0.9;  ///< Lambda reduction parameter for conservative backtracking

    floatingpoint BACKTRACKSLOPE = 0.4;   ///< Backtracking slope

    const floatingpoint QUADTOL = 0.1; // Used in Quadratic line search
    const floatingpoint ETOTALTOL = 1e-7;//Relative energy change tolerance.
    const floatingpoint LAMBDAQUADTOL = 1e-3; //if values change less than 0.1% between
    // successive runs, consider it converged.

    //additional parameters to help store additional parameters
    floatingpoint minimumE = (floatingpoint) 1e10;
    floatingpoint TotalEnergy = (floatingpoint)0.0;
    floatingpoint maxForcebackup = (floatingpoint)0.0;
    //@}

    // Track the past 100 lambdas.
    //@{
    uint maxprevlambdacount = 10;
    vector<floatingpoint> previouslambdavec=vector<floatingpoint>(maxprevlambdacount,0.0);
    short headpos = 0; //position where the next lambda can be inserted.
    short count = 0;//counter to track the number of successful lambda attempts made.
    float sum = 0;//sum of the lambdas found in previouslambdavcec.
    bool runningaveragestatus = false;
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


    double allFDotF();
	double allFADotFA();
	double allFADotFAP();
	double allFDotFA();
    
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
    void shiftGradient(double d);

    void setgradients(){
        FADotFA = allFADotFA();
	    FADotFAP = allFADotFAP();
    }
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
                                         floatingpoint maxForce,
                                         floatingpoint LAMBDAMAX,
                                         floatingpoint LAMBDARUNNINGAVERAGEPROBABILITY,
                                         bool *gpu_safestate, bool *ETolstate);
    
    /// The safemode backtracking search, returns the first energy decrease
    ///@note - The most robust linesearch method, but very slow

    floatingpoint safeBacktrackingLineSearch(
        ForceFieldManager& FFM, floatingpoint MAXDIST, floatingpoint maxForce,
        floatingpoint LAMBDAMAX, bool *gpu_safestate, bool *ETolstate);

    ///Quadratic line search introduced from LAMMPS based on Dennis and Schnabel
    floatingpoint quadraticLineSearch(ForceFieldManager& FFM, floatingpoint MAXDIST,
                                         floatingpoint maxForce,
                                         floatingpoint LAMBDAMAX,
                                         floatingpoint LAMBDARUNNINGAVERAGEPROBABILITY,
                                         bool *gpu_safestate, bool *ETolstate);

	floatingpoint safeBacktrackingLineSearchV2(ForceFieldManager& FFM, floatingpoint
										 MAXDIST, floatingpoint maxForce, floatingpoint
										 LAMBDAMAX, bool *gpu_safestate, bool *M_ETolstate);

	floatingpoint quadraticLineSearchV2(ForceFieldManager& FFM, floatingpoint MAXDIST,
	                                  floatingpoint maxForce,
	                                  floatingpoint LAMBDAMAX,
	                                  floatingpoint LAMBDARUNNINGAVERAGEPROBABILITY,
	                                  bool *gpu_safestate, bool *ETolstate);

	floatingpoint quadraticoptimization(ForceFieldManager& FFM, const
	vector<floatingpoint>& lambdavec,
			const vector<floatingpoint>&  energyvec);

    void setLAMBDATOL(int maxF_order){

        int orderdimension = 3; ///1000s of nm
        //Float gives you a minimum of 9 sig figs. If you operate in a 10^3nm system, the
        // decimal part can go upto 10^-6. In our system, Lambda*F_i should not be
        // greater than 10^-6. We would like to be cautious and ensure that all numbers
        // have  to the order of 10^-3. Hence, O(Lambda*F_i) >= 10^-(9-3-3) = 10^-3
        int LAMBDATOLorder = -(9-orderdimension -3) - maxF_order;
        LAMBDATOL = 1;
        if(LAMBDATOLorder > 0){
            for(int i =0; i < LAMBDATOLorder; i ++)
                LAMBDATOL *= 10;
        }
        else{
            for(int i =0; i > LAMBDATOLorder; i --)
                LAMBDATOL *= 0.1;
        }
		//Since, wthe force threshold are in 10^0 of pN at the lowest range, our lambda
		// should be in the order of 10^-5.
        LAMBDATOL = max<floatingpoint>(1e-8, LAMBDATOL);
        LAMBDATOL = min<floatingpoint>(1e-5, LAMBDATOL);

//        cout<<"maxF order "<<maxF_order<<" lambdatol "<<LAMBDATOL<<endl;
    }

    void setBACKTRACKSLOPE(floatingpoint _btslope){BACKTRACKSLOPE = _btslope;}

    void copycoordsifminimumE(floatingpoint maxForce){

    	if(TotalEnergy <= minimumE){
    		//update minimum energy
    		minimumE = TotalEnergy;
    		maxForcebackup = maxForce;
    		//take backup of coordinates.
		    const std::size_t num = Bead::getDbData().coords.size_raw();
		    Bead::getDbData().coords_minE.resize(num);
		    for(size_t i = 0; i < num; ++i) {
			    Bead::getDbData().coords_minE.value[i] = Bead::getDbData().coords.value[i];
		    }
    	}
    }

    void copybackupcoordinates(){

    	if(Bead::getDbData().coords_minE.size()) {
		    cout<<"Copying coordinates with the lowest energy during minimization "<<endl;
		    cout<<"Energy = "<<minimumE<<" pN.nm"<<endl;
		    cout<<"MaxForce = "<<maxForcebackup<<" pN "<<endl;
		    const std::size_t num = Bead::getDbData().coords.size_raw();
		    for (size_t i = 0; i < num; ++i) {
			    Bead::getDbData().coords.value[i] = Bead::getDbData().coords_minE.value[i];
		    }
	    }
    }

    //@}

    /// Print forces on all beads
    void printForces();
    
public:
    [[deprecated]] static long N; ///< Number of beads in the system, set before each minimization
    
    virtual ~CGMethod() {};

    /// Minimize the system
    virtual MinimizationResult minimize(ForceFieldManager &FFM, floatingpoint GRADTOL,
                          floatingpoint MAXDIST, floatingpoint LAMBDAMAX,
                          floatingpoint LAMBDARUNNINGAVERAGEPROBABILITY,
                          string _LINESEARCHALGORITHM,
                          bool steplimit) = 0;

};


#endif
