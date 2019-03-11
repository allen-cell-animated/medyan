

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

#include "CGPolakRibiereMethod.h"

#include "ForceFieldManager.h"
#include "Composite.h"
#include "Output.h"
#include "Bubble.h"
#include "cross_check.h"
#ifdef CUDAACCL
#include "nvToolsExt.h"
#endif
void PolakRibiere::minimize(ForceFieldManager &FFM, double GRADTOL,
                            double MAXDIST, double LAMBDAMAX, bool steplimit){
#ifdef CUDATIMETRACK
    chrono::high_resolution_clock::time_point tbeginTot, tendTot;
    chrono::high_resolution_clock::time_point tbeginII, tendII;
    tbeginTot = chrono::high_resolution_clock::now();
    tbeginII = chrono::high_resolution_clock::now();
#endif
    //number of steps
    int N;
    if(steplimit) {

        int beadMaxStep = 5 * Bead::numBeads();
        N = (beadMaxStep > _MINNUMSTEPS ? beadMaxStep : _MINNUMSTEPS);
    }
    else
        N = numeric_limits<int>::max();
    startMinimization();//TODO needs to be hostallocdefault and MemCpyAsync followed by CudaStreamSynchronize

    FFM.vectorizeAllForceFields();//each forcefield needs to use hostallocdefault and MemCpyAsync followed by CudaStreamSynchronize

    FFM.computeForces(coord, force); //split and synchronize in the end

#ifdef SERIAL // SERIAL
    FFM.copyForces(forceAux, force);
    FFM.copyForces(forceAuxPrev, force);
#endif
    
    //M as the first letter in variables signifies that it is used by minimizer
    // (as opposed to finding lambda)
    bool Ms_isminimizationstate, Ms_issafestate;
    int numIter = 0;
    double lambda;
#ifdef SERIAL
    //TODO Comment during SERIAL_CUDACROSSCHECK @{
    bool *Mc_isminimizationstate;
    bool *Mc_issafestate;
//    @}

    Mc_isminimizationstate = new bool[1];
    Mc_issafestate = new bool[1];
    Mc_isminimizationstate[0] = false;//points to address of Mmh_stop
    Mc_issafestate[0] = false;//points to address of Msh_stop

    //FIND MAXIMUM ERROR BETWEEN CUDA AND VECTORIZED FORCES{
    //VECTORIZED. Prep for Polak{
    double curGrad = CGMethod::allFDotF();
    Ms_isminimizationstate = true;
    Ms_issafestate = false;
    Ms_isminimizationstate = maxF() > GRADTOL;
    //}
    //
#ifdef DETAILEDOUTPUT
    std::cout<<"printing beads & forces"<<endl;
    long i = 0;
    long index = 0;
    for(auto b:Bead::getBeads()){
        index = 3 * b->_dbIndex;
        std::cout<<b->getID()<<" "<<coord[index]<<" "<<coord[index + 1]<<" "
                ""<<coord[index + 2]<<" "
                ""<<force[index]<<" "
                ""<<force[index + 1]<<" "<<force[index + 2]<<endl;
    }
    std::cout<<"printed beads & forces"<<endl;
#endif
    //
#endif
    
    numIter = 0;

//std::cout<<"----------------------------------------"<<endl;
//    std::cout<<"maxF "<<maxF()<<endl;
#ifdef SERIAL
    while (/* Iteration criterion */  numIter < N &&
           /* Gradient tolerance  */  (Ms_isminimizationstate )) {
//#ifdef CUDATIMETRACK_MACRO
//        chrono::high_resolution_clock::time_point tbeginiter, tenditer;
//        tbeginiter = chrono::high_resolution_clock::now();
//#endif

        double beta, newGrad, prevGrad;
//        std::cout<<"SERL maxF "<<maxF()<<endl;

        numIter++;


        bool *dummy;
        lambda = _safeMode ? safeBacktrackingLineSearch(FFM, MAXDIST, LAMBDAMAX, dummy)
                           : backtrackingLineSearch(FFM, MAXDIST, LAMBDAMAX, dummy);

        if(Ms_isminimizationstate)
            //SERIAL VERSION
            moveBeads(lambda);

        //compute new forces
        FFM.computeForces(coord, forceAux);//split and synchronize
#ifdef DETAILEDOUTPUT
        std::cout<<"MB printing beads & forces L "<<lambda<<endl;
        long i = 0;
        long index = 0;
        for(auto b:Bead::getBeads()){
            index = 3 * b->_dbIndex;

            std::cout<<b->getID()<<" "<<coord[index]<<" "<<coord[index + 1]<<" "
                    ""<<coord[index + 2]<<" "
                    ""<<forceAux[index]<<" "
                    ""<<forceAux[index + 1]<<" "<<forceAux[index + 2]<<" "<<3 *
                    b->_dbIndex<<endl;
        }
        std::cout<<"MB printed beads & forces"<<endl;
#endif

        //compute direction
//        std::cout<<"serial"<<endl;
        newGrad = CGMethod::allFADotFA();
        prevGrad = CGMethod::allFADotFAP();


        //Polak-Ribieri update
        beta = max(0.0, (newGrad - prevGrad) / curGrad);
        if(Ms_isminimizationstate)
            //shift gradient
            shiftGradient(beta);

        //vectorized copy
        FFM.copyForces(forceAuxPrev, forceAux);

        //direction reset if not downhill, or no progress made
        Ms_issafestate = CGMethod::allFDotFA() <= 0 || areEqual(curGrad, newGrad);
        if(Ms_issafestate && Ms_isminimizationstate ) {
            shiftGradient(0.0);
            _safeMode = true;

        }

        curGrad = newGrad;
        auto maxForce = maxF();
        Ms_isminimizationstate = maxForce > GRADTOL;
        
        for(auto b : Bubble::getBubbles()) {
            cout << "test_end = " << b->getBead()->coordinate[2] << endl;
        }

    }
#endif //SERIAL
//    std::cout<<"SERL Total number of iterations "<<numIter<<endl;

    if (numIter >= N) {
#ifdef CUDAACCL
//        FFM.CUDAcopyForces(*Msp1, CUDAcommon::getCUDAvars().gpu_forceAux,CUDAcommon::getCUDAvars().gpu_force);
//        //copy back forces and calculate load forces in CPU.
//        CUDAcommon::handleerror(cudaDeviceSynchronize());
#endif
        //TODO think about integrating CUDA version here
        cout << endl;

        cout << "WARNING: Did not minimize in N = " << N << " steps." << endl;
        cout << "Maximum force in system = " << maxF() << endl;

        cout << "Culprit ..." << endl;
        auto b = maxBead();
        if(b != nullptr) b->getParent()->printSelf();

#ifdef CUDAACCL
        auto cvars = CUDAcommon::getCUDAvars();
        cvars.streamvec.clear();
        CUDAcommon::cudavars = cvars;
#endif
        cout << "System energy..." << endl;
        FFM.computeEnergy(coord, force, 0.0, true);
#ifdef CUDAACCL
        for(auto strm:CUDAcommon::getCUDAvars().streamvec)
            CUDAcommon::handleerror(cudaStreamSynchronize(*strm));
#endif
        cout << endl;
    }

    //final force calculation
    FFM.computeForces(coord, force);

#ifdef SERIAL
    FFM.copyForces(forceAux, force);
#endif

#ifdef SERIAL
    //TODO Comment during SERIAL_CUDACROSSCHECK @{
    delete [] Mc_isminimizationstate;
    delete [] Mc_issafestate;
    //@}
#endif
 
    endMinimization();
    FFM.computeLoadForces();
    //std::cout<<"End Minimization************"<<endl;
    FFM.cleanupAllForceFields();

#ifdef DETAILEDOUTPUT
    std::cout<<"printing beads & forces"<<endl;
    for(auto b:Bead::getBeads()){
        std::cout<<b->getID()<<" "<<b->coordinate[0]<<" "<<b->coordinate[1]<<" "
                ""<<b->coordinate[2]<<" "
                ""<<b->force[index]<<" "
                ""<<b->force[index+1]<<" "<<b->force[index+2]<<endl;
    }
    std::cout<<"printed beads & forces"<<endl;
#endif

}
