

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
#ifdef CUDAACCL
#include "nvToolsExt.h"
#endif
void PolakRibiere::minimize(ForceFieldManager &FFM, double GRADTOL,
                            double MAXDIST, double LAMBDAMAX, bool steplimit){

    //number of steps
    int N;
    if(steplimit) {
        int beadMaxStep = 5 * Bead::numBeads();
        N = (beadMaxStep > _MINNUMSTEPS ? beadMaxStep : _MINNUMSTEPS);
//        N = 100;
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

#ifdef SERIAL // SERIAL
    FFM.copyForces(forceAux, force);
    FFM.copyForces(forceAuxPrev, force);
#endif
    //M as the first letter in variables signifies that it is used by minimizer
    // (as opposed to finding lambda)
    bool Ms_isminimizationstate, Ms_issafestate;
    int numIter = 0;
    double lambda;
#ifdef CUDAACCL
    volatile bool *Mc_isminimizationstate;
    volatile bool *Mc_issafestate;
    Ms_isminimizationstate = false;
    Ms_issafestate = false;
#endif
#ifdef SERIAL
    //TODO Comment during SERIAL_CUDACROSSCHECK @{
    bool *Mc_isminimizationstate;
    bool *Mc_issafestate;
    //@}

    Mc_isminimizationstate = new bool[1];
    Mc_issafestate = new bool[1];
    Mc_isminimizationstate[0] = false;//points to address of Mmh_stop
    Mc_issafestate[0] = false;//points to address of Msh_stop
#endif

#ifdef  CUDAACCL
    //wait for forces to be calculated
//    nvtxRangePushA("CCFsync");
    for(auto strm:CUDAcommon::getCUDAvars().streamvec)
        CUDAcommon::handleerror(cudaStreamSynchronize(*strm));
//    nvtxRangePop();

//    nvtxRangePushA("CCPFstream");
    if(!(CUDAcommon::getCUDAvars().conservestreams) || stream1 == NULL)
        CUDAcommon::handleerror(cudaStreamCreate(&stream1));
    if(!(CUDAcommon::getCUDAvars().conservestreams) || stream1 == NULL)
        CUDAcommon::handleerror(cudaStreamCreate(&stream2));
    if(!(CUDAcommon::getCUDAvars().conservestreams) || stream1 == NULL)
        CUDAcommon::handleerror(cudaStreamCreate(&stream3));
//    nvtxRangePop();
//    nvtxRangePushA("CCPF");
    FFM.CUDAcopyForces(stream1, CUDAcommon::getCUDAvars().gpu_forceAux,CUDAcommon::getCUDAvars().gpu_force);//pass a
    // stream
    FFM.CUDAcopyForces(stream2, CUDAcommon::getCUDAvars().gpu_forceAuxP,CUDAcommon::getCUDAvars().gpu_force);//pass a
    // stream
//    nvtxRangePop();
    double *gpu_GRADTOL;
    double gradtol[1];
    gradtol[0]= GRADTOL;
    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_GRADTOL, sizeof(double)));
    CUDAcommon::handleerror(cudaMemcpy(gpu_GRADTOL, gradtol, sizeof(double), cudaMemcpyHostToDevice));
    CGMethod::CUDAallFDotF(stream3);//curGrad //pass a stream
    //synchronize streams
//    nvtxRangePushA("CCPFstream");
    CUDAcommon::handleerror(cudaStreamSynchronize(stream1));
    CUDAcommon::handleerror(cudaStreamSynchronize(stream2));
    CUDAcommon::handleerror(cudaStreamSynchronize(stream3));
    if(!(CUDAcommon::getCUDAvars().conservestreams)) {
        CUDAcommon::handleerror(cudaStreamDestroy(stream1));
        CUDAcommon::handleerror(cudaStreamDestroy(stream2));
        CUDAcommon::handleerror(cudaStreamDestroy(stream3));
    }
//    nvtxRangePop();
//PING PONG
    bool  *Mmh_stop, *Mmg_stop1, *Mmg_stop2, *Mmg_s1, *Mmg_s2, *Mmg_ss;//minimization state
    bool  *Msh_stop, *Msg_stop1, *Msg_stop2, *Msg_s1, *Msg_s2, *Msg_ss;//safe state

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

//    nvtxRangePushA("MEvcreate");
    if(!(CUDAcommon::getCUDAvars().conservestreams) || Ms1 == NULL)
        CUDAcommon::handleerror(cudaStreamCreate(&Ms1));
    if(!(CUDAcommon::getCUDAvars().conservestreams) || Ms2 == NULL)
        CUDAcommon::handleerror(cudaStreamCreate(&Ms2));
    if(!(CUDAcommon::getCUDAvars().conservestreams) || Ms3 == NULL)
        CUDAcommon::handleerror(cudaStreamCreate(&Ms3));
    if(!(CUDAcommon::getCUDAvars().conservestreams) || Ms4 == NULL)
        CUDAcommon::handleerror(cudaStreamCreate(&Ms4));
    if(!(CUDAcommon::getCUDAvars().conservestreams) || Me1 == NULL)
        CUDAcommon::handleerror(cudaEventCreate(&Me1));
    if(!(CUDAcommon::getCUDAvars().conservestreams) || Me2 == NULL)
        CUDAcommon::handleerror(cudaEventCreate(&Me2));
    if(!(CUDAcommon::getCUDAvars().conservestreams) || stream_shiftsafe == NULL)
        CUDAcommon::handleerror(cudaStreamCreate(&stream_shiftsafe));
    if(!(CUDAcommon::getCUDAvars().conservestreams) || stream_dotcopy == NULL)
        CUDAcommon::handleerror(cudaStreamCreate(&stream_dotcopy));
    if(!(CUDAcommon::getCUDAvars().conservestreams) || event_safe == NULL)
        CUDAcommon::handleerror(cudaEventCreate(&event_safe));
    if(!(CUDAcommon::getCUDAvars().conservestreams) || event_dot == NULL)
        CUDAcommon::handleerror(cudaEventCreate(&event_dot));
//    nvtxRangePop();

    Mmh_stop[0] = true; //Minimizationstate //Yes = Minimize. No = Don't minimize.
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

//    nvtxRangePushA("PolakCudacalc");
//    CGMethod::CUDAgetPolakvars(false, *Msp1, gpu_GRADTOL, Mmg_s1, Mmg_s2, Msg_s2, Mc_isminimizationstate);
    std::cout<<endl;
    CGMethod::CUDAgetPolakvars(*Msp1, gpu_GRADTOL, Mmg_s1, Mmg_s2, Mc_isminimizationstate);
    CUDAcommon::handleerror(cudaEventRecord(*Mep1, *Msp1));
    CUDAcommon::handleerror(cudaGetLastError(),"CUDAgetPolakvars", "CGPolakRibiereMethod.cu");
#ifdef SERIAL_CUDACROSSCHECK
    CUDAcommon::handleerror(cudaDeviceSynchronize());
    std::cout<<"FMAX SL "<<maxF()<<endl;
#endif
//    nvtxRangePop();
    //Copy to host
//    nvtxRangePushA("PolakCudawait");
    CUDAcommon::handleerror(cudaStreamWaitEvent(Ms3, *Mep1, 0));
//    nvtxRangePop();
//    nvtxRangePushA("PolakCudacopy");
    cudaMemcpyAsync(Mmh_stop, Mmg_s2, sizeof(bool), cudaMemcpyDeviceToHost, Ms3);
//    nvtxRangePop();
//@}
#endif
#ifdef SERIAL //SERIAL
    //FIND MAXIMUM ERROR BETWEEN CUDA AND VECTORIZED FORCES{
    //VECTORIZED. Prep for Polak{
    double curGrad = CGMethod::allFDotF();
    Ms_isminimizationstate = true;
    Ms_issafestate = false;
//    nvtxRangePushA("Polakserial");
    Ms_isminimizationstate = maxF() > GRADTOL;
//    nvtxRangePop();
    //}
#endif

    while (/* Iteration criterion */  numIter < N &&
           /* Gradient tolerance  */  (Ms_isminimizationstate ||
           Mc_isminimizationstate[0])) {
#ifdef CUDAACCL
//        //@{
//        size_t free1, total;
//        CUDAcommon::handleerror(cudaMemGetInfo(&free1, &total));
//        cudaFree(0);
//        CUDAcommon::handleerror(cudaMemGetInfo(&free1, &total));
//        std::cout<<"Polak Before iter "<<numIter<<". Free VRAM in bytes "<<free1<<". Total VRAM in bytes "
//                 <<total<<endl;
//        auto cvars = CUDAcommon::getCUDAvars();
//        cvars.memincuda  = free1;
//        CUDAcommon::cudavars = cvars;
//        //@}
//PING PONG SWAP
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
#endif
#ifdef SERIAL_CUDACROSSCHECK
        std::cout<<"******************"<<endl;
#endif
#ifdef SERIAL
        double beta, newGrad, prevGrad;
        std::cout<<"SL maxF "<<maxF()<<endl;
#endif
//PING ENDS
        numIter++;
#ifdef SERIAL_CUDACROSSCHECK
        std::cout<<"SL safestate "<<_safeMode<<endl;
#endif
#ifdef CUDAACCL
        if(Mc_issafestate[0]) {
            _safeMode = false;
        }
//        nvtxRangePushA("while_Polak_sync");
//        CUDAcommon::handleerror(cudaStreamSynchronize(*Msp2));//make sure previous iteration is done.
//        nvtxRangePop();
        //find lambda by line search, move beads
//        nvtxRangePushA("lambda");
        lambda = backtrackingLineSearchCUDA(FFM, MAXDIST, LAMBDAMAX, Msg_s1);
//        nvtxRangePop();
#endif

#ifdef SERIAL
        bool *dummy;
//        nvtxRangePushA("lambda");
        lambda = _safeMode ? safeBacktrackingLineSearch(FFM, MAXDIST, LAMBDAMAX, dummy)
                           : backtrackingLineSearch(FFM, MAXDIST, LAMBDAMAX, dummy);
//        nvtxRangePop();
#endif
#ifdef SERIAL_CUDACROSSCHECK

        CUDAcommon::handleerror(cudaDeviceSynchronize());
        double cuda_lambda[1];
        CUDAcommon::handleerror(cudaMemcpy(cuda_lambda, CUDAcommon::cudavars.gpu_lambda,  sizeof(double),
                                           cudaMemcpyDeviceToHost));
        std::cout<<"Lambda "<<cuda_lambda[0]<<" "<<lambda<<endl;
#endif
#ifdef CUDAACCL
//        std::cout<<"move beads"<<endl;
        CUDAmoveBeads(*Msp1, Mmg_s1);
//        CUDAcommon::handleerror(cudaEventRecord(*Mep1, *Msp1));//This seems unnecessary.
        //wait for movebeads to finish before calculating forces1
        CUDAcommon::handleerror(cudaStreamSynchronize(*Msp1));
#endif
#ifdef SERIAL
        if(Ms_isminimizationstate)
            //SERIAL VERSION
            moveBeads(lambda);
#endif
#if defined(CROSSCHECK) || defined(CUDAACCL)
        cross_checkclass::Aux=true;
#endif
#ifdef CUDAACCL
        cvars = CUDAcommon::getCUDAvars();
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
//        nvtxRangePushA("CCPF"); //Copy forces
        FFM.CUDAcopyForces(stream_dotcopy, CUDAcommon::getCUDAvars().gpu_forceAuxP,CUDAcommon::getCUDAvars().gpu_forceAux);
//        nvtxRangePop();
        //Polak-Ribieri update beta & shift gradient
//        std::cout<<"shift Gradient"<<endl;
        CUDAshiftGradient(stream_shiftsafe, Mmg_s1);
#ifdef SERIAL_CUDACROSSCHECK
        CUDAcommon::handleerror(cudaStreamSynchronize(stream_shiftsafe));
#endif
        //@CUDA Get minimizaton state{
//        nvtxRangePushA("Polakcudacalc");
//        CGMethod::CUDAgetPolakvars(true, *Msp1, gpu_GRADTOL, Mmg_s1, Mmg_s2, Msg_s2, Mc_isminimizationstate);

        CGMethod::CUDAgetPolakvars(*Msp1, gpu_GRADTOL, Mmg_s1, Mmg_s2, Mc_isminimizationstate);
        CUDAcommon::handleerror(cudaEventRecord(*Mep1, *Msp1));
        CGMethod::CUDAgetPolakvars2(stream_shiftsafe, Msg_s2);

//        nvtxRangePop();
//        std::cout<<"shift Gradient Safe"<<endl;
#ifdef SERIAL_CUDACROSSCHECK
        CUDAcommon::handleerror(cudaDeviceSynchronize());
        std::cout<<"FMAX SL "<<maxF()<<endl;
#endif
        CUDAshiftGradientifSafe(stream_shiftsafe, Mmg_s1, Msg_s1);
//        CUDAcommon::handleerror(cudaDeviceSynchronize());

        if(Mc_isminimizationstate[0]  == true){
            //Copy to host
//            nvtxRangePushA("Polakcudawait");
            CUDAcommon::handleerror(cudaStreamWaitEvent(Ms3, *Mep1, 0));
//            CUDAcommon::handleerror(cudaStreamWaitEvent(Ms4, event_safe, 0));//event_safe is not attached to any event

//            nvtxRangePop();
//            nvtxRangePushA("Polakcudacopy");
//            std::cout<<"min state copy"<<endl;
            cudaMemcpyAsync(Mmh_stop, Mmg_s2, sizeof(bool), cudaMemcpyDeviceToHost, Ms3);
            //TODO remove later... June 7, 2018. Removed.
//            std::cout<<"safe state copy"<<endl;
//            cudaMemcpyAsync(Msh_stop, Msg_s2, sizeof(bool), cudaMemcpyDeviceToHost, Ms4);
//            nvtxRangePop();
        }
#endif
#ifdef SERIAL
        //compute direction
//        std::cout<<"serial"<<endl;
        newGrad = CGMethod::allFADotFA();
        prevGrad = CGMethod::allFADotFAP();

        //Polak-Ribieri update
        beta = max(0.0, (newGrad - prevGrad) / curGrad);
        if(Ms_isminimizationstate)
            //shift gradient
            shiftGradient(beta);
#ifdef SERIAL_CUDACROSSCHECK
        CUDAcommon::handleerror(cudaDeviceSynchronize(),"CGPolakRibiereMethod.cu","CGPolakRibiereMethod.cu");
        std::cout<<"Beta serial "<<beta<<endl;
        std::cout<<"newGrad "<<newGrad<<" prevGrad "<<prevGrad<<" curGrad "<<curGrad<<endl;
#endif
        //vectorized copy
//        nvtxRangePushA("SCPF");
//        std::cout<<"copy forces serial"<<endl;
        FFM.copyForces(forceAuxPrev, forceAux);
//        nvtxRangePop();
        //direction reset if not downhill, or no progress made
//        nvtxRangePushA("Polakserial");
        Ms_issafestate = CGMethod::allFDotFA() <= 0 || areEqual(curGrad, newGrad);
//        nvtxRangePop();
        if(Ms_issafestate && Ms_isminimizationstate ) {
            shiftGradient(0.0);
            _safeMode = true;
        }
        curGrad = newGrad;
//        nvtxRangePushA("Polakserial");
        Ms_isminimizationstate = maxF() > GRADTOL;
//        nvtxRangePop();
#endif

//        std::cout<<"M "<<Mc_isminimizationstate[0]<<" "<<Ms_isminimizationstate<<endl;
//        std::cout<<endl;
    }
    std::cout<<"Total iterations "<<numIter<<endl;
//    std::cout<<"maxF "<<maxF()<<" "<<GRADTOL<<" "<<Ms_isminimizationstate<<" "<<Mc_isminimizationstate[0]<<endl;

    if (numIter >= N) {
#ifdef CUDAACCL
//        FFM.CUDAcopyForces(*Msp1, CUDAcommon::getCUDAvars().gpu_forceAux,CUDAcommon::getCUDAvars().gpu_force);
//        //copy back forces and calculate load forces in CPU.
//        CUDAcommon::handleerror(cudaDeviceSynchronize());
#endif
        //TODO think about integrating CUDA version here
        cout << endl;

        cout << "WARNING: Did not minimize in N = " << N << " steps." << endl;
//        cout << "Maximum force in system = " << maxF() << endl;

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
//        nvtxRangePushA("CCFsync");
        for(auto strm:CUDAcommon::getCUDAvars().streamvec)
            CUDAcommon::handleerror(cudaStreamSynchronize(*strm));
//        nvtxRangePop();
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
//    nvtxRangePushA("CCFsync");
    for(auto strm:CUDAcommon::getCUDAvars().streamvec)
        CUDAcommon::handleerror(cudaStreamSynchronize(*strm));
//    nvtxRangePop();
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
    if(!(CUDAcommon::getCUDAvars().conservestreams)) {
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
    }
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
#endif
#ifdef SERIAL
    //TODO Comment during SERIAL_CUDACROSSCHECK @{
    delete [] Mc_isminimizationstate;
    delete [] Mc_issafestate;
    //@}
#endif
    //TODO make sure it calculates stretchforce in CUDA.
#ifdef CUDAACCL
    FFM.assignallforcemags();
#endif
    endMinimization();
    FFM.computeLoadForces();
    std::cout<<"End minimization-----------------"<<endl;

    FFM.cleanupAllForceFields();
}
