

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
#include "cross_check.h"
#include "nvToolsExt.h"
void PolakRibiere::minimize(ForceFieldManager &FFM, double GRADTOL,
                            double MAXDIST, double LAMBDAMAX, bool steplimit){

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

#ifdef CUDAACCL
    cross_checkclass::Aux=false;
    auto cvars = CUDAcommon::getCUDAvars();
    cvars.streamvec.clear();
    CUDAcommon::cudavars = cvars;
#endif

    FFM.computeForces(coord, force); //split and synchronize in the end

#ifndef CUDAACCL // SERIAL
    nvtxRangePushA("SCPF");
    //
    FFM.copyForces(forceAux, force);
    FFM.copyForces(forceAuxPrev, force);
    //
    nvtxRangePop();
#endif
    //M as the first letter in variables signifies that it is used by minimizer (as opposed to finding lambda)
    bool Ms_isminimizationstate, Ms_issafestate;
    int numIter = 0;
    double lambda;
#ifdef CUDAACCL
    volatile bool *Mc_isminimizationstate;
    volatile bool *Mc_issafestate;
    Ms_isminimizationstate = false;
    Ms_issafestate = false;
#else
    bool *Mc_isminimizationstate;
    bool *Mc_issafestate;
    Mc_isminimizationstate = new bool[1];
    Mc_issafestate = new bool[1];
    Mc_isminimizationstate[0] = false;//points to address of Mmh_stop
    Mc_issafestate[0] = false;//points to address of Msh_stop
#endif

#ifdef  CUDAACCL
    //wait for forces to be calculated
    nvtxRangePushA("CCFsync");
    for(auto strm:CUDAcommon::getCUDAvars().streamvec)
        CUDAcommon::handleerror(cudaStreamSynchronize(*strm));
    nvtxRangePop();

    nvtxRangePushA("CCPFstream");
    cudaStream_t stream1, stream2, stream3;
    CUDAcommon::handleerror(cudaStreamCreate(&stream1));
    CUDAcommon::handleerror(cudaStreamCreate(&stream2));
    CUDAcommon::handleerror(cudaStreamCreate(&stream3));
    nvtxRangePop();
    nvtxRangePushA("CCPF");
    FFM.CUDAcopyForces(stream1, CUDAcommon::getCUDAvars().gpu_forceAux,CUDAcommon::getCUDAvars().gpu_force);//pass a
    // stream
    FFM.CUDAcopyForces(stream2, CUDAcommon::getCUDAvars().gpu_forceAuxP,CUDAcommon::getCUDAvars().gpu_force);//pass a
    // stream
    nvtxRangePop();
    double *gpu_GRADTOL;
    double gradtol[1];
    gradtol[0]= GRADTOL;
    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_GRADTOL, sizeof(double)));
    CUDAcommon::handleerror(cudaMemcpy(gpu_GRADTOL, gradtol, sizeof(double), cudaMemcpyHostToDevice));
    CGMethod::CUDAallFDotF(stream3);//curGrad //pass a stream
    //synchronize streams
    nvtxRangePushA("CCPFstream");
    CUDAcommon::handleerror(cudaStreamSynchronize(stream1));
    CUDAcommon::handleerror(cudaStreamSynchronize(stream2));
    CUDAcommon::handleerror(cudaStreamSynchronize(stream3));
    CUDAcommon::handleerror(cudaStreamDestroy(stream1));
    CUDAcommon::handleerror(cudaStreamDestroy(stream2));
    CUDAcommon::handleerror(cudaStreamDestroy(stream3));
    nvtxRangePop();
//PING PONG
    bool  *Mmh_stop, *Mmg_stop1, *Mmg_stop2, *Mmg_s1, *Mmg_s2, *Mmg_ss;//minimization state
    bool  *Msh_stop, *Msg_stop1, *Msg_stop2, *Msg_s1, *Msg_s2, *Msg_ss;//safe state
    cudaStream_t Ms1, Ms2, Ms3, Ms4, *Msp1, *Msp2, *Msps;
    cudaEvent_t Me1, Me2, *Mep1, *Mep2, *Meps;

    //PING PONG
    //minimization state
    cudaMalloc(&Mmg_stop1, sizeof(bool));
    cudaMalloc(&Mmg_stop2, sizeof(bool));
    cudaHostAlloc(&Mmh_stop, sizeof(bool), cudaHostAllocDefault);
    //safe state
    cudaMalloc(&Msg_stop1, sizeof(bool));
    cudaMalloc(&Msg_stop2, sizeof(bool));
    cudaHostAlloc(&Msh_stop, sizeof(bool), cudaHostAllocDefault);
    //@
    //Memory alloted
    //@{
//    size_t allocmem = 0;
//    allocmem += sizeof(double) + 4 * sizeof(bool);
//    auto c = CUDAcommon::getCUDAvars();
//    c.memincuda += allocmem;
//    CUDAcommon::cudavars = c;
//    std::cout<<"Total allocated memory "<<c.memincuda/1024<<endl;
//    std::cout<<"Memory allocated "<< allocmem/1024<<"Memory freed 0"<<endl;
    //@}
    nvtxRangePushA("MEvcreate");
    CUDAcommon::handleerror(cudaStreamCreate(&Ms1));
    CUDAcommon::handleerror(cudaStreamCreate(&Ms2));
    CUDAcommon::handleerror(cudaStreamCreate(&Ms3));
    CUDAcommon::handleerror(cudaStreamCreate(&Ms4));
    CUDAcommon::handleerror(cudaEventCreate(&Me1));
    CUDAcommon::handleerror(cudaEventCreate(&Me2));
    CUDAcommon::handleerror(cudaStreamCreate(&stream_shiftsafe));
    CUDAcommon::handleerror(cudaStreamCreate(&stream_dotcopy));
    CUDAcommon::handleerror(cudaEventCreate(&event_safe));
    CUDAcommon::handleerror(cudaEventCreate(&event_dot));
    nvtxRangePop();

    Mmh_stop[0] = true; //Minimizationstate
    Msh_stop[0] = false; //safe state
    Mc_isminimizationstate = Mmh_stop;//points to address of Mmh_stop
    Mc_issafestate = Msh_stop;//points to address of Msh_stop

    Msp1 = &Ms1;
    Msp2 = &Ms2;
    Mep1 = &Me1;
    Mep2 = &Me2;
    Mmg_s1 = Mmg_stop1;
    Mmg_s2 = Mmg_stop2;
    Msg_s1 = Msg_stop1;
    Msg_s2 = Msg_stop2;
//set Mmg_stop1, Mmg_stop2 to true and Msg_stop1, Msg_stop2 to false.

// stick to single stream going forward.
//@CUDA Get minimizaton state{
    //calculate MAXF
    CGMethod::CUDAinitializePolak(*Msp1, Mmg_s1, Mmg_s2, Msg_s1, Msg_s2);
    CUDAcommon::handleerror(cudaGetLastError(),"CUDAinitializePolak", "CGPolakRibiereMethod.cu");

    nvtxRangePushA("PolakCudacalc");
//    CGMethod::CUDAgetPolakvars(false, *Msp1, gpu_GRADTOL, Mmg_s1, Mmg_s2, Msg_s2, Mc_isminimizationstate);
    CGMethod::CUDAgetPolakvars(*Msp1, gpu_GRADTOL, Mmg_s1, Mmg_s2, Mc_isminimizationstate);
    CUDAcommon::handleerror(cudaEventRecord(*Mep1, *Msp1));
    CUDAcommon::handleerror(cudaGetLastError(),"CUDAgetPolakvars", "CGPolakRibiereMethod.cu");
    nvtxRangePop();
    //Copy to host
    nvtxRangePushA("PolakCudawait");
    CUDAcommon::handleerror(cudaStreamWaitEvent(Ms3, *Mep1, 0));
    nvtxRangePop();
    nvtxRangePushA("PolakCudacopy");
    cudaMemcpyAsync(Mmh_stop, Mmg_s2, sizeof(bool), cudaMemcpyDeviceToHost, Ms3);
    nvtxRangePop();
//@}
#else //SERIAL
    //FIND MAXIMUM ERROR BETWEEN CUDA AND VECTORIZED FORCES{
    //VECTORIZED. Prep for Polak{
    double curGrad = CGMethod::allFDotF();
    Ms_isminimizationstate = true;
    Ms_issafestate = false;
    nvtxRangePushA("Polakserial");
    Ms_isminimizationstate = maxF() > GRADTOL;
    nvtxRangePop();
    //}
#endif
    while (/* Iteration criterion */  numIter < N &&
                                      /* Gradient tolerance  */  (Ms_isminimizationstate ||
                                                                  Mc_isminimizationstate[0])) {
//PING PONG SWAP
#ifdef CUDAACCL
//        CUDAcommon::handleerror(cudaStreamWaitEvent(*Msp2, *Mep1, 0));
        CUDAcommon::handleerror(cudaStreamSynchronize(*Msp1));
        CUDAcommon::handleerror(cudaStreamSynchronize(stream_shiftsafe));
        Msps = Msp1;
        Msp1 = Msp2;
        Msp2 = Msps;
        Meps = Mep1;
        Mep1 = Mep2;
        Mep2 = Meps;
        Mmg_ss = Mmg_s1;
        Mmg_s1 = Mmg_s2;
        Mmg_s2 = Mmg_ss;
        Msg_ss = Msg_s1;
        Msg_s1 = Msg_s2;
        Msg_s2 = Msg_ss;
#else
        double beta, newGrad, prevGrad;
        std::cout<<"maxF "<<maxF()<<endl;
#endif
//PING ENDS
        numIter++;
#ifdef CUDAACCL
        if(Mc_issafestate[0]) {
            _safeMode = false;
        }
        nvtxRangePushA("while_Polak_sync");
//        CUDAcommon::handleerror(cudaStreamSynchronize(*Msp2));//make sure previous iteration is done.
        nvtxRangePop();
        //find lambda by line search, move beads
        nvtxRangePushA("lambda");
        lambda = _safeMode ? safeBacktrackingLineSearch(FFM, MAXDIST, LAMBDAMAX, Msg_s1)
                           : backtrackingLineSearch(FFM, MAXDIST, LAMBDAMAX, Msg_s1);
        nvtxRangePop();
#else
        bool *dummy;
        nvtxRangePushA("lambda");
        lambda = _safeMode ? safeBacktrackingLineSearch(FFM, MAXDIST, LAMBDAMAX, dummy)
                           : backtrackingLineSearch(FFM, MAXDIST, LAMBDAMAX, dummy);
        nvtxRangePop();
#endif

#ifdef CUDAACCL
//        std::cout<<"move beads"<<endl;
        CUDAmoveBeads(*Msp1, Mmg_s1);
        CUDAcommon::handleerror(cudaEventRecord(*Mep1, *Msp1));
        CUDAcommon::handleerror(cudaGetLastError(),"CUDAmoveBeads", "CGPolakRibiereMethod.cu");
        //wait for movebeads to finish before calculating forces
        CUDAcommon::handleerror(cudaStreamSynchronize(*Msp1));
#else
        if(Ms_isminimizationstate)
            //SERIAL VERSION
            moveBeads(lambda);
#endif

#if defined(CROSSCHECK) || defined(CUDAACCL)
        cross_checkclass::Aux=true;
#endif
#ifdef CUDAACCL
        auto cvars = CUDAcommon::getCUDAvars();
        cvars.streamvec.clear();
        CUDAcommon::cudavars = cvars;
        CUDAcommon::handleerror(cudaStreamSynchronize(stream_dotcopy));
#endif
        //compute new forces
//        std::cout<<"compute forces"<<endl;
        FFM.computeForces(coord, forceAux);//split and synchronize
#ifdef  CUDAACCL
        //wait for forces to be calculated
        for(auto strm:CUDAcommon::getCUDAvars().streamvec)
            CUDAcommon::handleerror(cudaStreamSynchronize(*strm));

        //compute direction CUDA
//        std::cout<<"FdotFA"<<endl;
        CGMethod::CUDAallFADotFA(stream_dotcopy); //newGrad
//        CUDAcommon::handleerror(cudaDeviceSynchronize());
//        std::cout<<"FdotFAP"<<endl;
        CGMethod::CUDAallFADotFAP(stream_dotcopy); //prevGrad
//        CUDAcommon::handleerror(cudaDeviceSynchronize());
//        std::cout<<"copy forces"<<endl;
        CUDAcommon::handleerror(cudaEventRecord(event_dot,stream_dotcopy));
        CUDAcommon::handleerror(cudaStreamWaitEvent(stream_shiftsafe, event_dot,0));
        nvtxRangePushA("CCPF"); //Copy forces
        FFM.CUDAcopyForces(stream_dotcopy, CUDAcommon::getCUDAvars().gpu_forceAuxP,CUDAcommon::getCUDAvars().gpu_forceAux);
        nvtxRangePop();
        //Polak-Ribieri update beta & shift gradient
//        std::cout<<"shift Gradient"<<endl;
        CUDAshiftGradient(stream_shiftsafe, Mmg_s1);
        //@CUDA Get minimizaton state{
        nvtxRangePushA("Polakcudacalc");
//        CGMethod::CUDAgetPolakvars(true, *Msp1, gpu_GRADTOL, Mmg_s1, Mmg_s2, Msg_s2, Mc_isminimizationstate);
        CGMethod::CUDAgetPolakvars(*Msp1, gpu_GRADTOL, Mmg_s1, Mmg_s2, Mc_isminimizationstate);
        CUDAcommon::handleerror(cudaEventRecord(*Mep1, *Msp1));
        CGMethod::CUDAgetPolakvars2(stream_shiftsafe, Msg_s2);
        CUDAcommon::handleerror(cudaGetLastError(),"CUDAgetPolakvars", "CGPolakRibiereMethod.cu");
        nvtxRangePop();
//        std::cout<<"shift Gradient Safe"<<endl;
        CUDAshiftGradientifSafe(stream_shiftsafe, Mmg_s1, Msg_s1);
        CUDAcommon::handleerror(cudaGetLastError(),"CUDAshiftGradientifSafe", "CGPolakRibiereMethod.cu");
//        CUDAcommon::handleerror(cudaDeviceSynchronize());
        if(Mc_isminimizationstate[0]  == true){
            //Copy to host
            nvtxRangePushA("Polakcudawait");
            CUDAcommon::handleerror(cudaStreamWaitEvent(Ms3, *Mep1, 0));
            CUDAcommon::handleerror(cudaStreamWaitEvent(Ms4, event_safe, 0));
            nvtxRangePop();
            nvtxRangePushA("Polakcudacopy");
//            std::cout<<"min state copy"<<endl;
            cudaMemcpyAsync(Mmh_stop, Mmg_s2, sizeof(bool), cudaMemcpyDeviceToHost, Ms3);
            //TODO remove later.
//            std::cout<<"safe state copy"<<endl;
            cudaMemcpyAsync(Msh_stop, Msg_s2, sizeof(bool), cudaMemcpyDeviceToHost, Ms4);
            nvtxRangePop();
        }

#else
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
        nvtxRangePushA("SCPF");
//        std::cout<<"copy forces serial"<<endl;
        FFM.copyForces(forceAuxPrev, forceAux);
        nvtxRangePop();
        //direction reset if not downhill, or no progress made
        nvtxRangePushA("Polakserial");
        Ms_issafestate = CGMethod::allFDotFA() <= 0 || areEqual(curGrad, newGrad);
        nvtxRangePop();
        if(Ms_issafestate && Ms_isminimizationstate ) {
            shiftGradient(0.0);
            _safeMode = true;
        }
        curGrad = newGrad;
        nvtxRangePushA("Polakserial");
        Ms_isminimizationstate = maxF() > GRADTOL;
        nvtxRangePop();
#endif
//        std::cout<<"M "<<Mc_isminimizationstate[0]<<" "<<Ms_isminimizationstate<<endl;
//        std::cout<<endl;
    }
    std::cout<<"Total iterations "<<numIter<<endl;
//    std::cout<<"maxF "<<maxF()<<" "<<GRADTOL<<" "<<Ms_isminimizationstate<<" "<<Mc_isminimizationstate[0]<<endl;

    if (numIter >= N) {
#ifdef CUDAACCL
        CUDAcommon::handleerror(cudaDeviceSynchronize());
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
        nvtxRangePushA("CCFsync");
        for(auto strm:CUDAcommon::getCUDAvars().streamvec)
            CUDAcommon::handleerror(cudaStreamSynchronize(*strm));
        nvtxRangePop();
#endif
        cout << endl;
    }

//    cout << "Minimized." << endl;

#if defined(CROSSCHECK) || defined(CUDAACCL)
    cross_checkclass::Aux=false;
#endif
#ifdef CUDAACCL
    cvars = CUDAcommon::getCUDAvars();
    cvars.streamvec.clear();
    CUDAcommon::cudavars = cvars;
#endif

    //final force calculation
    FFM.computeForces(coord, force);
#ifdef CUDAACCL
    nvtxRangePushA("CCFsync");
    for(auto strm:CUDAcommon::getCUDAvars().streamvec)
        CUDAcommon::handleerror(cudaStreamSynchronize(*strm));
    nvtxRangePop();
#endif
#ifdef CUDAACCL
    FFM.CUDAcopyForces(*Msp1, CUDAcommon::getCUDAvars().gpu_forceAux,CUDAcommon::getCUDAvars().gpu_force);
    //copy back forces and calculate load forces in CPU.
#endif

    FFM.copyForces(forceAux, force);
#ifdef CUDAACCL
    CUDAcommon::handleerror(cudaFreeHost(Mmh_stop));
    CUDAcommon::handleerror(cudaFree(Mmg_stop1));
    CUDAcommon::handleerror(cudaFree(Mmg_stop2));
    CUDAcommon::handleerror(cudaFree(Msg_stop1));
    CUDAcommon::handleerror(cudaFree(Msg_stop2));
    CUDAcommon::handleerror(cudaFree(gpu_GRADTOL));
    CUDAcommon::handleerror(cudaFreeHost(Msh_stop));
    CUDAcommon::handleerror(cudaStreamSynchronize(Ms1));
    CUDAcommon::handleerror(cudaStreamSynchronize(Ms2));
    CUDAcommon::handleerror(cudaStreamSynchronize(Ms3));
    CUDAcommon::handleerror(cudaStreamSynchronize(Ms4));
    CUDAcommon::handleerror(cudaStreamDestroy(Ms1));
    CUDAcommon::handleerror(cudaStreamDestroy(Ms2));
    CUDAcommon::handleerror(cudaStreamDestroy(Ms3));
    CUDAcommon::handleerror(cudaStreamDestroy(Ms4));
    CUDAcommon::handleerror(cudaEventDestroy(Me1));
    CUDAcommon::handleerror(cudaEventDestroy(Me2));
    CUDAcommon::handleerror(cudaEventDestroy(event_safe));
    CUDAcommon::handleerror(cudaEventDestroy(event_dot));
    CUDAcommon::handleerror(cudaStreamDestroy(stream_dotcopy));
    CUDAcommon::handleerror(cudaStreamDestroy(stream_shiftsafe));
    //Memory alloted
    //@{
//    allocmem = 0;
//    allocmem += sizeof(double) + 4 * sizeof(bool);
//    c = CUDAcommon::getCUDAvars();
//    c.memincuda -= allocmem;
//    CUDAcommon::cudavars = c;
//    std::cout<<"Total allocated memory "<<c.memincuda/1024<<endl;
//    std::cout<<"Memory allocated 0 . Memory freed "<<allocmem/1024<<endl;
    //@}
#else
    delete Mc_isminimizationstate;
    delete Mc_issafestate;
#endif
    endMinimization();
    FFM.computeLoadForces();
    std::cout<<"End minimization-----------------"<<endl;

    FFM.cleanupAllForceFields();
}
