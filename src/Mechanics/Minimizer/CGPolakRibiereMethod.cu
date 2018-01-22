

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
    else {
        N = numeric_limits<int>::max();
    }

    startMinimization();//needs to be hostallocdefault and MemCpyAsync followed by CudaStreamSynchronize
    FFM.vectorizeAllForceFields();//each forcefield needs to use hostallocdefault and MemCpyAsync followed by CudaStreamSynchronize

#if defined(CROSSCHECK) || defined(CUDAENABLED) ||defined(CUDACROSSCHECK)
    cross_checkclass::Aux=false;
#endif

#ifdef CUDAACCL
    auto cvars = CUDAcommon::getCUDAvars();
    cvars.streamvec.clear();
    CUDAcommon::cudavars = cvars;
#endif

    FFM.computeForces(coord, force); //split and synchronize in the end

//    bool isminimizestate = true;
//    bool issafestate = false;
    nvtxRangePushA("SCPF");
    //
    FFM.copyForces(forceAux, force);
    FFM.copyForces(forceAuxPrev, force);
    //
    nvtxRangePop();

#ifdef  CUDAACCL
    //wait for forces to be calculated
    nvtxRangePushA("CCFsync");
    for(auto strm:CUDAcommon::getCUDAvars().streamvec)
        CUDAcommon::handleerror(cudaStreamSynchronize(*strm));
    nvtxRangePop();

//    bool* boolvars =  new bool[2];
//    //boolvars[0] = maxF() > GRADTOL
//    // boolvars[1] = CGMethod::allFDotFA() <= 0 || areEqual(curGrad, newGrad)
//    boolvars[0] = false;
//    boolvars[1] = false;
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

    //M as the first letter in variables signifies that it is used by minimizer (as opposed to finding lambda)
    volatile bool *Mc_isminimizationstate;
    volatile bool *Mc_issafestate;
    bool Ms_isminimizationstate, Ms_issafestate;
    bool *Mmh_stop, *Mmg_stop1, *Mmg_stop2, *Mmg_s1, *Mmg_s2, *Mmg_ss;//minimization state
    bool *Msh_stop, *Msg_stop1, *Msg_stop2, *Msg_s1, *Msg_s2, *Msg_ss;//safe state
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

    nvtxRangePushA("MEvcreate");
    CUDAcommon::handleerror(cudaStreamCreate(&Ms1));
    CUDAcommon::handleerror(cudaStreamCreate(&Ms2));
    CUDAcommon::handleerror(cudaStreamCreate(&Ms3));
    CUDAcommon::handleerror(cudaStreamCreate(&Ms4));
    CUDAcommon::handleerror(cudaEventCreate(&Me1));
    CUDAcommon::handleerror(cudaEventCreate(&Me2));
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
    CUDAcommon::handleerror(cudaGetLastError());

    nvtxRangePushA("PolakCudacalc");
    CGMethod::CUDAgetPolakvars(false, *Msp1, gpu_GRADTOL, Mmg_s1, Mmg_s2, Msg_s2, Mc_isminimizationstate);
    CUDAcommon::handleerror(cudaEventRecord(*Mep1, *Msp1));
    CUDAcommon::handleerror(cudaGetLastError());
    nvtxRangePop();
    //Copy to host
    nvtxRangePushA("PolakCudawait");
    CUDAcommon::handleerror(cudaStreamWaitEvent(Ms3, *Mep1, 0));
    nvtxRangePop();
    nvtxRangePushA("PolakCudacopy");
    cudaMemcpyAsync(Mmh_stop, Mmg_s2, sizeof(bool), cudaMemcpyDeviceToHost, Ms3);
    nvtxRangePop();
//@}
#endif

//FIND MAXIMUM ERROR BETWEEN CUDA AND VECTORIZED FORCES{
    double F_i[3*Bead::getBeads().size()];
    double max2;
    CUDAcommon::handleerror(cudaMemcpy(F_i, CUDAcommon::getCUDAvars().gpu_force, 3 * Bead::getBeads().size() *
                                                                                 sizeof(double), cudaMemcpyDeviceToHost));
    max2=0.0;
    for(int iter=0;iter<Bead::getBeads().size();iter++) {
//        std::cout << "C " << F_i[3 * iter] << " " << F_i[3 * iter + 1] << " " << F_i[3 * iter + 2] <<" ";
//        std::cout << "V "<<force[3 * iter] << " " << force[3 * iter + 1] << " " << force[3 * iter + 2] << endl;
        if(abs(F_i[3 * iter]-force[3 * iter])>max2)
            max2 = abs(F_i[3 * iter]-force[3 * iter]);
        else if(abs(F_i[3 * iter +1]-force[3 * iter +1])>max2)
            max2 = abs(F_i[3 * iter +1]-force[3 * iter +1]);
        else if(abs(F_i[3 * iter +2]-force[3 * iter +2])>max2)
            max2 = abs(F_i[3 * iter +1]-force[3 * iter +1]);
    }
    std::cout<<"force alone start "<<max2<<endl;
//}mAX Error
//Error    for(auto b:Bead::getBeads()) {
//        b->updatePosition();
//    }
    //VECTORIZED. Prep for Polak{
    double curGrad = CGMethod::allFDotF();
    Ms_isminimizationstate = true;
    Ms_issafestate = false;
    int numIter = 0;
    nvtxRangePushA("Polakserial");
    Ms_isminimizationstate = maxF() > GRADTOL;
    nvtxRangePop();
    //}
    std::cout<<maxF()<<" "<<GRADTOL<<endl;
    //Start Polak-Ribiere
    while (/* Iteration criterion */  numIter < N &&
                                      /* Gradient tolerance  */  (Ms_isminimizationstate ||
                                              Mc_isminimizationstate[0])) {
        std::cout<<maxF()<<" "<<GRADTOL<<endl;
//        for(auto b:Bead::getBeads()) {
//            b->updatePosition();
//        }
//PING PONG SWAP
#ifdef CUDAACCL
        CUDAcommon::handleerror(cudaStreamWaitEvent(*Msp2, *Mep1, 0));
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
//PING ENDS
        numIter++;
        double lambda, beta, newGrad, prevGrad;

#if defined(CROSSCHECK)
        auto state=cross_check::crosscheckforces(force);
#endif
//        nvtxRangePushA("Polakcudasafewait");
//        CUDAcommon::handleerror(cudaStreamSynchronize(Ms4));
//        nvtxRangePop();
        std::cout<<"S "<<Mc_issafestate[0]<<" "<<Ms_issafestate<<endl;
        if(Mc_issafestate[0] ) {
            _safeMode = true;
        }
        //find lambda by line search, move beads
#ifdef CUDAACCL
        lambda = _safeMode ? safeBacktrackingLineSearch(FFM, MAXDIST, LAMBDAMAX, Msg_s1)
                           : backtrackingLineSearch(FFM, MAXDIST, LAMBDAMAX, Msg_s1);
#else
        bool *dummy;
        lambda = _safeMode ? safeBacktrackingLineSearch(FFM, MAXDIST, LAMBDAMAX, dummy)
                           : backtrackingLineSearch(FFM, MAXDIST, LAMBDAMAX, dummy);
#endif
        //Copy gpu_lambda to see if they match.
        //TODO REMOVE LATER
        double cudalambda[1];
        CUDAcommon::handleerror(cudaMemcpy(cudalambda, CUDAcommon::getCUDAvars().gpu_lambda, sizeof(double),
                                           cudaMemcpyDeviceToHost));
        std::cout<<"lambda "<<lambda<<" "<<cudalambda[0]<<endl;

//        double *F_i, *C_i;
//        F_i = new double[3*Bead::getBeads().size()];
//        C_i = new double[3*Bead::getBeads().size()];
//
////        double F_i[3*Bead::getBeads().size()];
////        double C_i[3*Bead::getBeads().size()];
//        double max1,max2,max3;
//        CUDAcommon::handleerror(cudaMemcpy(C_i, CUDAcommon::getCUDAvars().gpu_coord, 3 * Bead::getBeads().size() *
//                                                                                     sizeof(double), cudaMemcpyDeviceToHost));
//
//        //Find Max error in COORDS betwen CUDA and vectorized versions{
//        max1=0.0;
//        for(int iter=0;iter<Bead::getBeads().size();iter++) {
////            std::cout << "C " << C_i[3 * iter] << " " << C_i[3 * iter + 1] << " " << C_i[3 * iter + 2] <<" ";
////            std::cout << "V "<<coord[3 * iter] << " " << coord[3 * iter + 1] << " " <<coord[3 * iter + 2] << endl;
//            if(abs(C_i[3 * iter]-coord[3 * iter])>max1)
//                max1 = abs(C_i[3 * iter]-coord[3 * iter]);
//            else if(abs(C_i[3 * iter +1]-coord[3 * iter +1])>max1)
//                max1 = abs(C_i[3 * iter +1]-coord[3 * iter +1]);
//            else if(abs(C_i[3 * iter +2]-coord[3 * iter +2])>max1)
//                max1 = abs(C_i[3 * iter +1]-coord[3 * iter +1]);
//        }
//
//        //}END


#ifdef CUDAACCL
        std::cout<<"move beads"<<endl;
        CUDAmoveBeads(*Msp1, Mmg_s1);
        CUDAcommon::handleerror(cudaEventRecord(*Mep1, *Msp1));
        CUDAcommon::handleerror(cudaGetLastError());
        //wait for movebeads to finish before calculating forces
        CUDAcommon::handleerror(cudaStreamSynchronize(*Msp1));
#endif
        if(Ms_isminimizationstate)
            //SERIAL VERSION
            moveBeads(lambda);

//Find Max error in COORDS & FORCES betwen CUDA and vectorized versions
//
//        max2=0.0;max3=0.0;
////        std::cout<<"forces"<<endl;
//        for(int iter=0;iter<Bead::getBeads().size();iter++) {
////            std::cout << "C " << F_i[3 * iter] << " " << F_i[3 * iter + 1] << " " << F_i[3 * iter + 2] <<" ";
////            std::cout << "V "<<force[3 * iter] << " " << force[3 * iter + 1] << " " << force[3 * iter + 2] << endl;
//            if(abs(F_i[3 * iter]-force[3 * iter])>max2)
//                max2 = abs(F_i[3 * iter]-force[3 * iter]);
//            else if(abs(F_i[3 * iter +1]-force[3 * iter +1])>max2)
//                max2 = abs(F_i[3 * iter +1]-force[3 * iter +1]);
//            else if(abs(F_i[3 * iter +2]-force[3 * iter +2])>max2)
//                max2 = abs(F_i[3 * iter +1]-force[3 * iter +1]);
//        }
////        std::cout<<"coord"<<endl;
//        for(int iter=0;iter<Bead::getBeads().size();iter++) {
////            std::cout << "C " << C_i[3 * iter] << " " << C_i[3 * iter + 1] << " " << C_i[3 * iter + 2] <<" ";
////            std::cout << "V "<<coord[3 * iter] << " " << coord[3 * iter + 1] << " " <<coord[3 * iter + 2] << endl;
//            if(abs(C_i[3 * iter]-coord[3 * iter])>max3)
//                max3 = abs(C_i[3 * iter]-coord[3 * iter]);
//            else if(abs(C_i[3 * iter +1]-coord[3 * iter +1])>max3)
//                max3 = abs(C_i[3 * iter +1]-coord[3 * iter +1]);
//            else if(abs(C_i[3 * iter +2]-coord[3 * iter +2])>max3)
//                max3 = abs(C_i[3 * iter +1]-coord[3 * iter +1]);
//        }
//        std::cout<<"check "<<max1<<" "<<max2<<" "<<max3<<endl;
////}END coord & forces error


#if defined(CROSSCHECK) || defined(CUDAENABLED) ||defined(CUDACROSSCHECK)
        cross_checkclass::Aux=true;
#endif
#ifdef CUDAACCL
        auto cvars = CUDAcommon::getCUDAvars();
        cvars.streamvec.clear();
        CUDAcommon::cudavars = cvars;
#endif
        //compute new forces
        std::cout<<"compute forces"<<endl;
        FFM.computeForces(coord, forceAux);//split and synchronize
//vectorized copy
        nvtxRangePushA("SCPF");
        std::cout<<"copy forces serial"<<endl;
        FFM.copyForces(forceAuxPrev, forceAux);
        nvtxRangePop();

#ifdef  CUDAACCL
        //wait for forces to be calculated
        std::cout<<"copy forces"<<endl;
        for(auto strm:CUDAcommon::getCUDAvars().streamvec)
            CUDAcommon::handleerror(cudaStreamSynchronize(*strm));
        nvtxRangePushA("CCPF");
        FFM.CUDAcopyForces(*Msp1, CUDAcommon::getCUDAvars().gpu_forceAuxP,CUDAcommon::getCUDAvars().gpu_forceAux);
        nvtxRangePop();
        //compute direction CUDA
        std::cout<<"FdotFA"<<endl;
        CGMethod::CUDAallFADotFA(*Msp1);
        std::cout<<"FdotFAP"<<endl;
        CGMethod::CUDAallFADotFAP(*Msp1);

        //Polak-Ribieri update beta & shift gradient
        std::cout<<"shift Gradient"<<endl;
        CUDAshiftGradient(*Msp1, Mmg_s1);

        //@CUDA Get minimizaton state{
        nvtxRangePushA("Polakcudacalc");
        CGMethod::CUDAgetPolakvars(true, *Msp1, gpu_GRADTOL, Mmg_s1, Mmg_s2, Msg_s2, Mc_isminimizationstate);
        CUDAcommon::handleerror(cudaEventRecord(*Mep1, *Msp1));
        CUDAcommon::handleerror(cudaGetLastError());
        nvtxRangePop();
        std::cout<<"shift Gradient Safe"<<endl;
        CUDAshiftGradientifSafe(*Msp1, Mmg_s1, Msg_s1);
        CUDAcommon::handleerror(cudaEventRecord(*Mep1, *Msp1));
        CUDAcommon::handleerror(cudaGetLastError());

        if(Mc_isminimizationstate[0]  == true){
            //Copy to host
            nvtxRangePushA("Polakcudawait");
            CUDAcommon::handleerror(cudaStreamWaitEvent(Ms3, *Mep1, 0));
            CUDAcommon::handleerror(cudaStreamWaitEvent(Ms4, *Mep1, 0));
            nvtxRangePop();
            nvtxRangePushA("Polakcudacopy");
            std::cout<<"min state copy"<<endl;
            cudaMemcpyAsync(Mmh_stop, Mmg_s2, sizeof(bool), cudaMemcpyDeviceToHost, Ms3);
            //TODO remove later.
            std::cout<<"safe state copy"<<endl;
            cudaMemcpyAsync(Msh_stop, Msg_s2, sizeof(bool), cudaMemcpyDeviceToHost, Ms4);
            nvtxRangePop();
        }

#endif
        //compute direction
        std::cout<<"serial"<<endl;
        newGrad = CGMethod::allFADotFA();
        prevGrad = CGMethod::allFADotFAP();

        //Polak-Ribieri update
        beta = max(0.0, (newGrad - prevGrad) / curGrad);
        if(Ms_isminimizationstate)
            //shift gradient
            shiftGradient(beta);
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
        std::cout<<"M "<<Mc_isminimizationstate[0]<<" "<<Ms_isminimizationstate<<endl;
    }

    if (numIter >= N) {
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

#if defined(CROSSCHECK) || defined(CUDAENABLED) || defined(CUDACROSSCHECK)
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
    CUDAcommon::handleerror(cudaMemcpy(F_i, CUDAcommon::getCUDAvars().gpu_force, 3 * Bead::getBeads().size() *
                                                                                 sizeof(double), cudaMemcpyDeviceToHost));
    max2=0.0;
    for(int iter=0;iter<Bead::getBeads().size();iter++) {
//        std::cout << F_i[3 * iter] << " " << F_i[3 * iter + 1] << " " << F_i[3 * iter + 2] <<" ";
//        std::cout <<force[3 * iter] << " " << force[3 * iter + 1] << " " << force[3 * iter + 2] << endl;
        if(abs(F_i[3 * iter]-force[3 * iter])>max2)
            max2 = abs(F_i[3 * iter]-force[3 * iter]);
        else if(abs(F_i[3 * iter +1]-force[3 * iter +1])>max2)
            max2 = abs(F_i[3 * iter +1]-force[3 * iter +1]);
        else if(abs(F_i[3 * iter +2]-force[3 * iter +2])>max2)
            max2 = abs(F_i[3 * iter +1]-force[3 * iter +1]);
    }

    std::cout<<"force alone end "<<max2<<endl;

#if defined(CROSSCHECK)
        auto state=cross_check::crosscheckforces(force);
#endif
#ifdef CUDAACCL
    FFM.CUDAcopyForces(*Msp1, CUDAcommon::getCUDAvars().gpu_forceAux,CUDAcommon::getCUDAvars().gpu_force);
    //copy back forces and calculate load forces in CPU.
#endif

    FFM.copyForces(forceAux, force);

    endMinimization();
    FFM.computeLoadForces();
std::cout<<"-----------------"<<endl;

    FFM.cleanupAllForceFields();
}
