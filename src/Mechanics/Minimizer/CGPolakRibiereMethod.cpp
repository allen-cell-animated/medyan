

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

#include "CGPolakRibiereMethod.h"

#include "ForceFieldManager.h"
#include "Composite.h"
#include "Output.h"
#include "cross_check.h"
#ifdef CUDAACCL
#include "nvToolsExt.h"
#endif
#include "Structure/Bead.h"

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
#ifdef ALLSYNC
    cudaDeviceSynchronize();
#endif
    FFM.vectorizeAllForceFields();//each forcefield needs to use hostallocdefault and MemCpyAsync followed by CudaStreamSynchronize
#ifdef ALLSYNC
    cudaDeviceSynchronize();
#endif
#ifdef CUDAACCL
    cross_checkclass::Aux=false;
    auto cvars = CUDAcommon::getCUDAvars();
    cvars.streamvec.clear();
    CUDAcommon::cudavars = cvars;
#endif
#ifdef ALLSYNC
    cudaDeviceSynchronize();
#endif
    FFM.computeForces(Bead::getDbData().coords.data(), Bead::getDbData().forces.data()); //split and synchronize in the end

#ifdef SERIAL // SERIAL
    Bead::getDbData().forcesAux = Bead::getDbData().forces;
    Bead::getDbData().forcesAuxP = Bead::getDbData().forces;

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
//    @}

    Mc_isminimizationstate = new bool[1];
    Mc_issafestate = new bool[1];
    Mc_isminimizationstate[0] = false;//points to address of Mmh_stop
    Mc_issafestate[0] = false;//points to address of Msh_stop
#endif
#ifdef CUDATIMETRACK
    chrono::high_resolution_clock::time_point tbegin, tend;
#endif
#ifdef  CUDAACCL
#ifdef CUDATIMETRACK
//    chrono::high_resolution_clock::time_point tbegin, tend;
    tbegin = chrono::high_resolution_clock::now();
#endif
#ifdef ALLSYNC
    cudaDeviceSynchronize();
#endif
    for(auto strm:CUDAcommon::getCUDAvars().streamvec)
        CUDAcommon::handleerror(cudaStreamSynchronize(*strm));

#ifdef CUDATIMETRACK
    tend= chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed_run(tend - tbegin);
    CUDAcommon::cudatime.TveccomputeF.push_back(elapsed_run.count());
    CUDAcommon::cudatime.TcomputeF += elapsed_run.count();
#endif
#ifdef CUDATIMETRACK

    //Reset lambda time tracker.
    CUDAcommon::cudatime.Tlambda = 0.0;
#endif

    if(!(CUDAcommon::getCUDAvars().conservestreams) || stream1 == NULL)
        CUDAcommon::handleerror(cudaStreamCreate(&stream1));
    if(!(CUDAcommon::getCUDAvars().conservestreams) || stream1 == NULL)
        CUDAcommon::handleerror(cudaStreamCreate(&stream2));
    if(!(CUDAcommon::getCUDAvars().conservestreams) || stream1 == NULL)
        CUDAcommon::handleerror(cudaStreamCreate(&stream3));
    FFM.CUDAcopyForces(stream1, CUDAcommon::getCUDAvars().gpu_forceAux,CUDAcommon::getCUDAvars().gpu_force);//pass a
    // stream
#ifdef ALLSYNC
    cudaDeviceSynchronize();
#endif
    FFM.CUDAcopyForces(stream2, CUDAcommon::getCUDAvars().gpu_forceAuxP,CUDAcommon::getCUDAvars().gpu_force);//pass a
    // stream
#ifdef ALLSYNC
    cudaDeviceSynchronize();
#endif
    double *gpu_GRADTOL;
    double gradtol[1];
    gradtol[0]= GRADTOL;
    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_GRADTOL, sizeof(double)));
    CUDAcommon::handleerror(cudaMemcpy(gpu_GRADTOL, gradtol, sizeof(double), cudaMemcpyHostToDevice));
    CGMethod::CUDAallFDotF(stream3);//curGrad //pass a stream

    //synchronize streams
    CUDAcommon::handleerror(cudaStreamSynchronize(stream1));
    CUDAcommon::handleerror(cudaStreamSynchronize(stream2));
    CUDAcommon::handleerror(cudaStreamSynchronize(stream3));
    if(!(CUDAcommon::getCUDAvars().conservestreams)) {
        CUDAcommon::handleerror(cudaStreamDestroy(stream1));
        CUDAcommon::handleerror(cudaStreamDestroy(stream2));
        CUDAcommon::handleerror(cudaStreamDestroy(stream3));
    }
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
    //calculate MAXF and allFDotFA
    CGMethod::CUDAinitializePolak(*Msp1, Mmg_s1, Mmg_s2, Msg_s1, Msg_s2);
#ifdef ALLSYNC
    cudaDeviceSynchronize();
#endif
    CGMethod::CUDAgetPolakvars(*Msp1, gpu_GRADTOL, Mmg_s1, Mmg_s2, Mc_isminimizationstate);
#ifdef ALLSYNC
    cudaDeviceSynchronize();
#endif
    CUDAcommon::handleerror(cudaEventRecord(*Mep1, *Msp1));
    CUDAcommon::handleerror(cudaGetLastError(),"CUDAgetPolakvars", "CGPolakRibiereMethod.cu");
#ifdef SERIAL_CUDACROSSCHECK
    CUDAcommon::handleerror(cudaDeviceSynchronize());
    std::cout<<"FMAX SERL "<<maxF()<<endl;
#endif
#ifdef ALLSYNC
    cudaDeviceSynchronize();
#endif
    //Copy to host
    CUDAcommon::handleerror(cudaStreamWaitEvent(Ms3, *Mep1, 0));
    cudaMemcpyAsync(Mmh_stop, Mmg_s2, sizeof(bool), cudaMemcpyDeviceToHost, Ms3);
//@}
#endif
#ifdef SERIAL //SERIAL
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
        index = 3 * b->getIndex();
        std::cout<<b->getId()<<" "<< b->coordinate() <<" "
                "" << b->force() <<endl;
    }
    std::cout<<"printed beads & forces"<<endl;
#endif
    //
#endif
#ifdef CUDATIMETRACK_MACRO
    CUDAcommon::cudatime.Tlambda = 0.0;
    CUDAcommon::cudatime.Ecount = 0;
    CUDAcommon::cudatime.TcomputeE = 0.0;
    CUDAcommon::cudatime.TcomputeF = 0.0;
    double c1 = 0.0;
    double c2,c3,c4,c5;
    c2 = 0.0;c3 = 0.0;c4 = 0.0;c5=0.0;
#endif
#ifdef CUDATIMETRACK
    tendII= chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed_runslice1(tendII - tbeginII);
    std::cout<<"Slice time "<<elapsed_runslice1.count()<<endl;
    tbeginII = chrono::high_resolution_clock::now();
#endif
#ifdef CUDAACCL
    while (/* Iteration criterion */  numIter < N &&
           /* Gradient tolerance  */  (Mc_isminimizationstate[0])) {
//#ifdef CUDATIMETRACK_MACRO
//        chrono::high_resolution_clock::time_point tbeginiter, tenditer;
//        tbeginiter = chrono::high_resolution_clock::now();
//#endif
#ifdef CUDATIMETRACK
        chrono::high_resolution_clock::time_point tbegin, tend;
#endif

#ifdef ALLSYNC
        cudaDeviceSynchronize();
#endif
//PING PONG SWAP
//        CUDAcommon::handleerror(cudaStreamWaitEvent(*Msp2, *Mep1, 0));

#ifdef CUDATIMETRACK
        tbegin = chrono::high_resolution_clock::now();
#endif
        CUDAcommon::handleerror(cudaStreamSynchronize(*Msp1));
        CUDAcommon::handleerror(cudaStreamSynchronize(stream_shiftsafe));
#ifdef CUDATIMETRACK
        tend = chrono::high_resolution_clock::now();
        chrono::duration<double> elapsed_run1(tend - tbegin);
        c5 += elapsed_run1.count();
#endif

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
//PING ENDS
#ifdef ALLSYNC
        cudaDeviceSynchronize();
#endif
        numIter++;
#if defined(SERIAL_CUDACROSSCHECK) && defined(DETAILEDOUTPUT_LAMBDA)
        std::cout<<"SL safestate "<<_safeMode<<endl;
#endif
        if(Mc_issafestate[0]) {
            _safeMode = false;
        }
#ifdef CUDATIMETRACK_MACRO
        CUDAcommon::serltime.TcomputeEiter = 0.0;
        CUDAcommon::cudatime.TcomputeEiter = 0.0;
        chrono::high_resolution_clock::time_point tLbegin, tLend;
        tLbegin = chrono::high_resolution_clock::now();
#endif
//        CUDAcommon::handleerror(cudaStreamSynchronize(*Msp2));//make sure previous iteration is done.
        //find lambda by line search, move beads
        lambda = backtrackingLineSearchCUDA(FFM, MAXDIST, LAMBDAMAX, Msg_s1);
#ifdef CUDATIMETRACK_MACRO
        tLend= chrono::high_resolution_clock::now();
        chrono::duration<double> elapsed_runiter(tLend - tLbegin);
        std::cout<<" CUDA Iter "<<numIter<<" Lambda Time taken (s) "<<elapsed_runiter
                .count() - CUDAcommon::serltime.TcomputeEiter <<endl;
#endif
#ifdef ALLSYNC
        cudaDeviceSynchronize();
#endif

#ifdef CUDATIMETRACK
        tbegin = chrono::high_resolution_clock::now();
#endif
        CUDAmoveBeads(*Msp1, Mmg_s1);
//        CUDAcommon::handleerror(cudaEventRecord(*Mep1, *Msp1));//This seems unnecessary.
        //wait for movebeads to finish before calculating forces1

#ifdef ALLSYNC
        cudaDeviceSynchronize();
#endif
        //Synchronize to ensure forces are calculated on moved beads.
        CUDAcommon::handleerror(cudaStreamSynchronize(*Msp1));
//        CUDAcommon::handleerror(cudaStreamSynchronize(*Msp2));
        //@CUDA Get minimizaton state
        // @{
//        CGMethod::CUDAgetPolakvars(true, *Msp1, gpu_GRADTOL, Mmg_s1, Mmg_s2, Msg_s2, Mc_isminimizationstate);
        CGMethod::CUDAgetPolakvars(*Msp1, gpu_GRADTOL, Mmg_s1, Mmg_s2, Mc_isminimizationstate);
#ifdef ALLSYNC
        cudaDeviceSynchronize();
#endif
        CUDAcommon::handleerror(cudaEventRecord(*Mep1, *Msp1));
        //@}
#ifdef CUDATIMETRACK
        tend = chrono::high_resolution_clock::now();
        chrono::duration<double> elapsed_run1b(tend - tbegin);
        c1 += elapsed_run1b.count();
#endif

#if defined(CROSSCHECK) || defined(CUDAACCL)
        cross_checkclass::Aux=true;
#endif
#ifdef ALLSYNC
        cudaDeviceSynchronize();
#endif

        cvars = CUDAcommon::getCUDAvars();
        cvars.streamvec.clear();
        CUDAcommon::cudavars = cvars;
#ifdef CUDATIMETRACK
        tbegin = chrono::high_resolution_clock::now();
#endif
        CUDAcommon::handleerror(cudaStreamSynchronize(stream_dotcopy));
#ifdef CUDATIMETRACK
        tend = chrono::high_resolution_clock::now();
        chrono::duration<double> elapsed_run2(tend - tbegin);
        c2 += elapsed_run2.count();
#endif

        //compute new forces
        FFM.computeForces(Bead::getDbData().coords.data(), Bead::getDbData().forcesAux.data());//split and synchronize
#ifdef DETAILEDOUTPUT
        std::cout<<"MB printing beads & forces L "<<lambda<<endl;
        long i = 0;
        long index = 0;
        for(auto b:Bead::getBeads()){
            index = 3 * b->getIndex();

            std::cout<<b->getId()<<" "<<coord[index]<<" "<<coord[index + 1]<<" "
                    ""<<coord[index + 2]<<" "
                    ""<<forceAux[index]<<" "
                    ""<<forceAux[index + 1]<<" "<<forceAux[index + 2]<<" "<<3 *
                    b->getIndex()<<endl;
        }
        std::cout<<"MB printed beads & forces"<<endl;
#endif

#ifdef ALLSYNC
        cudaDeviceSynchronize();
#endif
#ifdef CUDATIMETRACK
        tbegin = chrono::high_resolution_clock::now();
#endif
        //wait for forces to be calculated
        for(auto strm:CUDAcommon::getCUDAvars().streamvec)
            CUDAcommon::handleerror(cudaStreamSynchronize(*strm));
#ifdef CUDATIMETRACK
        tend= chrono::high_resolution_clock::now();
        chrono::duration<double> elapsed_run(tend - tbegin);
        CUDAcommon::cudatime.TveccomputeF.push_back(elapsed_run.count());
        CUDAcommon::cudatime.TcomputeF += elapsed_run.count();
#endif
#ifdef CUDATIMETRACK
        tbegin = chrono::high_resolution_clock::now();
#endif
//Copying forces back to Host was here.
#ifdef CUDATIMETRACK
        tend = chrono::high_resolution_clock::now();
        chrono::duration<double> elapsed_run1c(tend - tbegin);
        c1 += elapsed_run1c.count();
#endif

#ifdef CUDATIMETRACK
        tbegin = chrono::high_resolution_clock::now();
#endif
        //compute direction CUDA
        CGMethod::CUDAallFADotFA(stream_dotcopy); //newGrad
#ifdef ALLSYNC
        cudaDeviceSynchronize();
#endif
        CGMethod::CUDAallFADotFAP(stream_dotcopy); //prevGrad
#ifdef ALLSYNC
        cudaDeviceSynchronize();
#endif
        CUDAcommon::handleerror(cudaEventRecord(event_dot,stream_dotcopy));
        CUDAcommon::handleerror(cudaStreamWaitEvent(stream_shiftsafe, event_dot,0));
//Copy forces
        FFM.CUDAcopyForces(stream_dotcopy, CUDAcommon::getCUDAvars().gpu_forceAuxP,CUDAcommon::getCUDAvars().gpu_forceAux);
#ifdef CUDATIMETRACK
        tend = chrono::high_resolution_clock::now();
        chrono::duration<double> elapsed_run2c(tend - tbegin);
        c2 += elapsed_run2c.count();
#endif
#ifdef ALLSYNC
        cudaDeviceSynchronize();
#endif
#ifdef CUDATIMETRACK
        tbegin = chrono::high_resolution_clock::now();
#endif
        //Polak-Ribieri update beta & shift gradient
        CUDAshiftGradient(stream_shiftsafe, Mmg_s1);
#ifdef SERIAL_CUDACROSSCHECK
        CUDAcommon::handleerror(cudaStreamSynchronize(stream_shiftsafe));
#endif
/*        CUDAcommon::handleerror(cudaStreamSynchronize(*Msp2));
        //@CUDA Get minimizaton state{
//        CGMethod::CUDAgetPolakvars(true, *Msp1, gpu_GRADTOL, Mmg_s1, Mmg_s2, Msg_s2, Mc_isminimizationstate);
        CGMethod::CUDAgetPolakvars(*Msp1, gpu_GRADTOL, Mmg_s1, Mmg_s2, Mc_isminimizationstate);
#ifdef ALLSYNC
        cudaDeviceSynchronize();
#endif
        CUDAcommon::handleerror(cudaEventRecord(*Mep1, *Msp1));*/
        CGMethod::CUDAgetPolakvars2(stream_shiftsafe, Msg_s2);
#ifdef ALLSYNC
        cudaDeviceSynchronize();
#endif
        CUDAshiftGradientifSafe(stream_shiftsafe, Mmg_s1, Msg_s1);
#ifdef ALLSYNC
        cudaDeviceSynchronize();
#endif
#ifdef CUDATIMETRACK
        tend = chrono::high_resolution_clock::now();
        chrono::duration<double> elapsed_run3(tend - tbegin);
        c3 += elapsed_run3.count();
#endif
#ifdef CUDATIMETRACK
        tbegin = chrono::high_resolution_clock::now();
#endif
        if(Mc_isminimizationstate[0]  == true){
            //Copy to host
            CUDAcommon::handleerror(cudaStreamWaitEvent(Ms3, *Mep1, 0));
//            CUDAcommon::handleerror(cudaStreamWaitEvent(Ms4, event_safe, 0));//event_safe is not attached to any event
#ifdef ALLSYNC
            cudaDeviceSynchronize();
#endif
            cudaMemcpyAsync(Mmh_stop, Mmg_s2, sizeof(bool), cudaMemcpyDeviceToHost, Ms3);
#ifdef ALLSYNC
            cudaDeviceSynchronize();
#endif
            //TODO remove later... June 7, 2018. Removed.
//            std::cout<<"safe state copy"<<endl;
//            cudaMemcpyAsync(Msh_stop, Msg_s2, sizeof(bool), cudaMemcpyDeviceToHost, Ms4);

        }
#ifdef CUDATIMETRACK
        tend = chrono::high_resolution_clock::now();
        chrono::duration<double> elapsed_run4(tend - tbegin);
        c4 += elapsed_run4.count();
#endif
    }

#ifdef CUDATIMETRACK_MACRO
    std::cout<<"CUDA start min time taken (s) "<<CUDAcommon::cudatime.Tstartmin<<endl;
    std::cout<<"CUDA Energy time taken (s) "<<CUDAcommon::cudatime.TcomputeE<<"for total "
            "iters "<<CUDAcommon::cudatime.Ecount<<endl;
    std::cout<<"CUDA Force time taken (s) "<<CUDAcommon::cudatime.TcomputeF<<endl;
    std::cout<<"CUDA Energy time per iter (s/iter) "<<CUDAcommon::cudatime.TcomputeE/
            (double(CUDAcommon::cudatime.Ecount))<<endl;
    std::cout<<"CUDA Split Times Iter "<<numIter<<" "<<c1<<" "<<c2<<" "<<c3<<" "<<c4<<" "
            ""<<c5<<endl;
    std::cout<<"CUDA Add "<<c5+c4+c3+c2+c1<<endl;
#endif
    std::cout<<"CUDA Total number of iterations "<<numIter<<endl;
#endif //CUDAACCL

    numIter = 0;
#ifdef CUDATIMETRACK_MACRO
    double s1 = 0.0;
    double s2,s3,s4;
    s2 = 0.0;s3 = 0.0;s4 = 0.0;
    CUDAcommon::serltime.Tlambda = 0.0;
    CUDAcommon::serltime.Ecount = 0;
    CUDAcommon::serltime.TcomputeE = 0.0;
    CUDAcommon::serltime.TcomputeF = 0.0;
#endif
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
#if defined(SERIAL_CUDACROSSCHECK) && defined(DETAILEDOUTPUT_LAMBDA)
        std::cout<<"SL safestate "<<_safeMode<<endl;
#endif
#ifdef CUDATIMETRACK_MACRO
        CUDAcommon::serltime.TcomputeEiter = 0.0;
        CUDAcommon::cudatime.TcomputeEiter = 0.0;
        chrono::high_resolution_clock::time_point tLbegin, tLend;
        tLbegin = chrono::high_resolution_clock::now();
#endif
        bool *dummy = nullptr;
        lambda = _safeMode ? safeBacktrackingLineSearch(FFM, MAXDIST, LAMBDAMAX, dummy)
                           : backtrackingLineSearch(FFM, MAXDIST, LAMBDAMAX, dummy);
#ifdef CUDATIMETRACK_MACRO
        tLend= chrono::high_resolution_clock::now();
        chrono::duration<double> elapsed_runiter(tLend - tLbegin);
        std::cout<<"SERL Iter "<<numIter<<" Lambda Time taken (s) "<<elapsed_runiter
                .count() - CUDAcommon::cudatime.TcomputeEiter
                 <<endl;
#endif

#ifdef SERIAL_CUDACROSSCHECK
        CUDAcommon::handleerror(cudaDeviceSynchronize());
        double cuda_lambda[1];
        CUDAcommon::handleerror(cudaMemcpy(cuda_lambda, CUDAcommon::cudavars.gpu_lambda,  sizeof(double),
                                           cudaMemcpyDeviceToHost));
        std::cout<<"Lambda CUDA "<<cuda_lambda[0]<<" SERL "<<lambda<<endl;
#endif
#ifdef CUDATIMETRACK
        tbegin = chrono::high_resolution_clock::now();
#endif
        if(Ms_isminimizationstate)
            //SERIAL VERSION
            moveBeads(lambda);

#if defined(CROSSCHECK) || defined(CUDAACCL)
        cross_checkclass::Aux=true;
#endif
#ifdef CUDATIMETRACK
        tend = chrono::high_resolution_clock::now();
        chrono::duration<double> elapsed_runs1(tend - tbegin);
        s1 += elapsed_runs1.count();
#endif
        //compute new forces
        FFM.computeForces(Bead::getDbData().coords.data(), Bead::getDbData().forcesAux.data());//split and synchronize
#ifdef DETAILEDOUTPUT
        std::cout<<"MB printing beads & forces L "<<lambda<<endl;
        long i = 0;
        long index = 0;
        for(auto b:Bead::getBeads()){
            index = 3 * b->getIndex();

            std::cout<<b->getId()<<" "<<coord[index]<<" "<<coord[index + 1]<<" "
                    ""<<coord[index + 2]<<" "
                    ""<<forceAux[index]<<" "
                    ""<<forceAux[index + 1]<<" "<<forceAux[index + 2]<<" "<<3 *
                    b->getIndex()<<endl;
        }
        std::cout<<"MB printed beads & forces"<<endl;
#endif
#ifdef CUDATIMETRACK
        tbegin = chrono::high_resolution_clock::now();
#endif
        //compute direction
//        std::cout<<"serial"<<endl;
        newGrad = CGMethod::allFADotFA();
        prevGrad = CGMethod::allFADotFAP();
#ifdef CUDATIMETRACK
        tend = chrono::high_resolution_clock::now();
        chrono::duration<double> elapsed_runs2a(tend - tbegin);
        s2 += elapsed_runs2a.count();
#endif
#ifdef CUDATIMETRACK
        tbegin = chrono::high_resolution_clock::now();
#endif
        //Polak-Ribieri update
        beta = max(0.0, (newGrad - prevGrad) / curGrad);
        if(Ms_isminimizationstate)
            //shift gradient
            shiftGradient(beta);
#ifdef CUDATIMETRACK
        tend = chrono::high_resolution_clock::now();
        chrono::duration<double> elapsed_runs3a(tend - tbegin);
        s3 += elapsed_runs3a.count();
#endif
#if defined(SERIAL_CUDACROSSCHECK) && defined(DETAILEDOUTPUT_BETA)
        std::cout<<"Shift Gradient "<<beta<<endl;
        CUDAcommon::handleerror(cudaDeviceSynchronize(),"CGPolakRibiereMethod.cu","CGPolakRibiereMethod.cu");
        std::cout<<"Beta serial "<<beta<<endl;
        std::cout<<"newGrad "<<newGrad<<" prevGrad "<<prevGrad<<" curGrad "<<curGrad<<endl;
#endif
#ifdef CUDATIMETRACK
        tbegin = chrono::high_resolution_clock::now();
#endif
#ifdef CUDATIMETRACK
        tbegin = chrono::high_resolution_clock::now();
#endif
        Bead::getDbData().forcesAuxP = Bead::getDbData().forcesAux;

#ifdef CUDATIMETRACK
        tend = chrono::high_resolution_clock::now();
        chrono::duration<double> elapsed_runs2b(tend - tbegin);
        s2 += elapsed_runs2b.count();
#endif
#ifdef CUDATIMETRACK
        tbegin = chrono::high_resolution_clock::now();
#endif
        //direction reset if not downhill, or no progress made
        Ms_issafestate = CGMethod::allFDotFA() <= 0 || areEqual(curGrad, newGrad);
        if(Ms_issafestate && Ms_isminimizationstate ) {
            shiftGradient(0.0);
            _safeMode = true;
#ifdef DETAILEDOUTPUT_LAMBDA
            std::cout<<"SERL FDOTFA "<<CGMethod::allFDotFA()<<" curGrad "<<curGrad<<" "
                    "newGrad "<<newGrad<<endl;
            std::cout<<"Shift Gradient 0.0"<<endl;
#endif
        }
#ifdef CUDATIMETRACK
        tend = chrono::high_resolution_clock::now();
        chrono::duration<double> elapsed_runs3b(tend - tbegin);
        s3 += elapsed_runs3b.count();
#endif
#ifdef CUDATIMETRACK
        tbegin = chrono::high_resolution_clock::now();
#endif
        curGrad = newGrad;
        auto maxForce = maxF();
        Ms_isminimizationstate = maxForce > GRADTOL;
#ifdef CUDATIMETRACK
        tend = chrono::high_resolution_clock::now();
        chrono::duration<double> elapsed_runs1b(tend - tbegin);
        s1 += elapsed_runs1b.count();
#endif
    }
#endif //SERIAL
//    std::cout<<"SERL Total number of iterations "<<numIter<<endl;
#ifdef CUDATIMETRACK_MACRO
    std::cout<<"SERL Energy time taken (s) "<<CUDAcommon::serltime.TcomputeE<<" for total "
            "iters "<<CUDAcommon::serltime.Ecount<<endl;
    std::cout<<"SERL Force time taken (s) "<<CUDAcommon::serltime.TcomputeF<<endl;
    std::cout<<"SERL Energy time per iter (s/iter) "<<CUDAcommon::serltime.TcomputeE/
                                                      (double(CUDAcommon::serltime.Ecount))
             <<endl;
    std::cout<<"SERL Split Times Iter "<<numIter<<" "<<s1<<" "<<s2<<" "<<s3<<" "<<s4<<endl;
    std::cout<<"SERL Add "<<s1+s2+s3+s4<<endl;
#endif
#ifdef CUDATIMETRACK
    tendII= chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed_runslice2(tendII - tbeginII);
    std::cout<<"Slice time "<<elapsed_runslice2.count()<<endl;
    tbeginII = chrono::high_resolution_clock::now();
#endif
//    while (/* Iteration criterion */  numIter < N &&
//           /* Gradient tolerance  */  (Ms_isminimizationstate ||
//           Mc_isminimizationstate[0])) {
//#ifdef CUDAACCL
//#ifdef ALLSYNC
//    cudaDeviceSynchronize();
//#endif
////PING PONG SWAP
////        CUDAcommon::handleerror(cudaStreamWaitEvent(*Msp2, *Mep1, 0));
//        CUDAcommon::handleerror(cudaStreamSynchronize(*Msp1));
//        CUDAcommon::handleerror(cudaStreamSynchronize(stream_shiftsafe));
//        Msps = Msp1;
//        Msp1 = Msp2;
//        Msp2 = Msps;
//        Meps = Mep1;
//        Mep1 = Mep2;
//        Mep2 = Meps;
//        Mmg_ss = Mmg_s1;
//        Mmg_s1 = Mmg_s2;
//        Mmg_s2 = Mmg_ss;
//        Msg_ss = Msg_s1;
//        Msg_s1 = Msg_s2;
//        Msg_s2 = Msg_ss;
//#endif
//#ifdef SERIAL
//        double beta, newGrad, prevGrad;
//        std::cout<<"SERL maxF "<<maxF()<<endl;
//#endif
////PING ENDS
//#ifdef ALLSYNC
//        cudaDeviceSynchronize();
//#endif
//        numIter++;
//#if defined(SERIAL_CUDACROSSCHECK) && defined(DETAILEDOUTPUT_LAMBDA)
//        std::cout<<"SL safestate "<<_safeMode<<endl;
//#endif
//#ifdef CUDAACCL
//        if(Mc_issafestate[0]) {
//            _safeMode = false;
//        }
////        CUDAcommon::handleerror(cudaStreamSynchronize(*Msp2));//make sure previous iteration is done.
//        //find lambda by line search, move beads
//        lambda = backtrackingLineSearchCUDA(FFM, MAXDIST, LAMBDAMAX, Msg_s1);
//#endif
//#ifdef ALLSYNC
//        cudaDeviceSynchronize();
//#endif
//#ifdef SERIAL
//        bool *dummy;
//        lambda = _safeMode ? safeBacktrackingLineSearch(FFM, MAXDIST, LAMBDAMAX, dummy)
//                           : backtrackingLineSearch(FFM, MAXDIST, LAMBDAMAX, dummy);
//#endif
//#ifdef SERIAL_CUDACROSSCHECK
//
//        CUDAcommon::handleerror(cudaDeviceSynchronize());
//        double cuda_lambda[1];
//        CUDAcommon::handleerror(cudaMemcpy(cuda_lambda, CUDAcommon::cudavars.gpu_lambda,  sizeof(double),
//                                           cudaMemcpyDeviceToHost));
//        std::cout<<"Lambda CUDA "<<cuda_lambda[0]<<" SERL "<<lambda<<endl;
//#endif
//#ifdef CUDAACCL
//        CUDAmoveBeads(*Msp1, Mmg_s1);
////        CUDAcommon::handleerror(cudaEventRecord(*Mep1, *Msp1));//This seems unnecessary.
//        //wait for movebeads to finish before calculating forces1
//#ifdef ALLSYNC
//        cudaDeviceSynchronize();
//#endif
//        CUDAcommon::handleerror(cudaStreamSynchronize(*Msp1));
//#endif
//#ifdef SERIAL
//        if(Ms_isminimizationstate)
//            //SERIAL VERSION
//            moveBeads(lambda);
//#endif
//#if defined(CROSSCHECK) || defined(CUDAACCL)
//        cross_checkclass::Aux=true;
//#endif
//#ifdef ALLSYNC
//        cudaDeviceSynchronize();
//#endif
//#ifdef CUDAACCL
//        cvars = CUDAcommon::getCUDAvars();
//        cvars.streamvec.clear();
//        CUDAcommon::cudavars = cvars;
//        CUDAcommon::handleerror(cudaStreamSynchronize(stream_dotcopy));
//#endif
//        //compute new forces
//        FFM.computeForces(coord, forceAux);//split and synchronize
//#ifdef DETAILEDOUTPUT
//        std::cout<<"MB printing beads & forces L "<<lambda<<endl;
//        long i = 0;
//        long index = 0;
//        for(auto b:Bead::getBeads()){
//            index = 3 * b->getIndex();
//
//            std::cout<<b->getId()<<" "<<coord[index]<<" "<<coord[index + 1]<<" "
//                    ""<<coord[index + 2]<<" "
//                    ""<<forceAux[index]<<" "
//                    ""<<forceAux[index + 1]<<" "<<forceAux[index + 2]<<" "<<3 *
//                    b->getIndex()<<endl;
//        }
//        std::cout<<"MB printed beads & forces"<<endl;
//#endif
//#ifdef  CUDAACCL
//#ifdef CUDATIMETRACK
//        chrono::high_resolution_clock::time_point tbegin, tend;
//        tbegin = chrono::high_resolution_clock::now();
//#endif
//#ifdef ALLSYNC
//        cudaDeviceSynchronize();
//#endif
//        //wait for forces to be calculated
//        for(auto strm:CUDAcommon::getCUDAvars().streamvec)
//            CUDAcommon::handleerror(cudaStreamSynchronize(*strm));
//#ifdef CUDATIMETRACK
//        tend= chrono::high_resolution_clock::now();
//        chrono::duration<double> elapsed_run(tend - tbegin);
//        CUDAcommon::cudatime.TveccomputeF.push_back(elapsed_run.count());
//        CUDAcommon::cudatime.TcomputeF += elapsed_run.count();
//#endif
//#ifdef CUDATIMETRACK
//        std::cout<<"Time total computeForces (s) CUDA "<<CUDAcommon::cudatime
//                .TcomputeF<<" SERL "<<CUDAcommon::serltime.TcomputeF<<" factor "
//                         ""<<CUDAcommon::serltime.TcomputeF/CUDAcommon::cudatime
//                .TcomputeF<<endl;
//        std::cout<<"Time split computeForces (s) CUDA ";
//        for(auto x:CUDAcommon::cudatime.TveccomputeF)
//            std::cout<<x<<" ";
//        std::cout<<endl;
//        std::cout<<"Time split computeForces (s) SERL ";
//        for(auto x:CUDAcommon::serltime.TveccomputeF)
//            std::cout<<x<<" ";
//        std::cout<<endl;
//#endif
//        //compute direction CUDA
//        CGMethod::CUDAallFADotFA(stream_dotcopy); //newGrad
//#ifdef ALLSYNC
//        cudaDeviceSynchronize();
//#endif
//        CGMethod::CUDAallFADotFAP(stream_dotcopy); //prevGrad
//#ifdef ALLSYNC
//        cudaDeviceSynchronize();
//#endif
//
//        CUDAcommon::handleerror(cudaEventRecord(event_dot,stream_dotcopy));
//        CUDAcommon::handleerror(cudaStreamWaitEvent(stream_shiftsafe, event_dot,0));
////Copy forces
//        FFM.CUDAcopyForces(stream_dotcopy, CUDAcommon::getCUDAvars().gpu_forceAuxP,CUDAcommon::getCUDAvars().gpu_forceAux);
//#ifdef ALLSYNC
//        cudaDeviceSynchronize();
//#endif
//        //Polak-Ribieri update beta & shift gradient
//        CUDAshiftGradient(stream_shiftsafe, Mmg_s1);
//#ifdef SERIAL_CUDACROSSCHECK
//        CUDAcommon::handleerror(cudaStreamSynchronize(stream_shiftsafe));
//#endif
//        //@CUDA Get minimizaton state{
////        CGMethod::CUDAgetPolakvars(true, *Msp1, gpu_GRADTOL, Mmg_s1, Mmg_s2, Msg_s2, Mc_isminimizationstate);
//        CGMethod::CUDAgetPolakvars(*Msp1, gpu_GRADTOL, Mmg_s1, Mmg_s2, Mc_isminimizationstate);
//#ifdef ALLSYNC
//        cudaDeviceSynchronize();
//#endif
//        CUDAcommon::handleerror(cudaEventRecord(*Mep1, *Msp1));
//        CGMethod::CUDAgetPolakvars2(stream_shiftsafe, Msg_s2);
//#ifdef ALLSYNC
//        cudaDeviceSynchronize();
//#endif
//#ifdef SERIAL_CUDACROSSCHECK
//        CUDAcommon::handleerror(cudaDeviceSynchronize());
//        std::cout<<"FMAX SERL "<<maxF()<<endl;
//#endif
//        CUDAshiftGradientifSafe(stream_shiftsafe, Mmg_s1, Msg_s1);
//#ifdef ALLSYNC
//        cudaDeviceSynchronize();
//#endif
//
//        if(Mc_isminimizationstate[0]  == true){
//            //Copy to host
//            CUDAcommon::handleerror(cudaStreamWaitEvent(Ms3, *Mep1, 0));
////            CUDAcommon::handleerror(cudaStreamWaitEvent(Ms4, event_safe, 0));//event_safe is not attached to any event
//#ifdef ALLSYNC
//            cudaDeviceSynchronize();
//#endif
//            cudaMemcpyAsync(Mmh_stop, Mmg_s2, sizeof(bool), cudaMemcpyDeviceToHost, Ms3);
//#ifdef ALLSYNC
//            cudaDeviceSynchronize();
//#endif
//            //TODO remove later... June 7, 2018. Removed.
////            std::cout<<"safe state copy"<<endl;
////            cudaMemcpyAsync(Msh_stop, Msg_s2, sizeof(bool), cudaMemcpyDeviceToHost, Ms4);
//        }
//#endif
//#ifdef SERIAL
//        //compute direction
////        std::cout<<"serial"<<endl;
//        newGrad = CGMethod::allFADotFA();
//        prevGrad = CGMethod::allFADotFAP();
//
//        //Polak-Ribieri update
////        std::cout<<newGrad<<" "<<prevGrad<<" "<<curGrad<<endl;
//        beta = max(0.0, (newGrad - prevGrad) / curGrad);
//        if(Ms_isminimizationstate)
//            //shift gradient
//            shiftGradient(beta);
//#if defined(SERIAL_CUDACROSSCHECK) && defined(DETAILEDOUTPUT_BETA)
//        std::cout<<"Shift Gradient "<<beta<<endl;
//        CUDAcommon::handleerror(cudaDeviceSynchronize(),"CGPolakRibiereMethod.cu","CGPolakRibiereMethod.cu");
//        std::cout<<"Beta serial "<<beta<<endl;
//        std::cout<<"newGrad "<<newGrad<<" prevGrad "<<prevGrad<<" curGrad "<<curGrad<<endl;
//#endif
//        //vectorized copy

////        std::cout<<"copy forces serial"<<endl;
//        FFM.copyForces(forceAuxPrev, forceAux);

//        //direction reset if not downhill, or no progress made

//        Ms_issafestate = CGMethod::allFDotFA() <= 0 || areEqual(curGrad, newGrad);

//        if(Ms_issafestate && Ms_isminimizationstate ) {
//            shiftGradient(0.0);
//            _safeMode = true;
//            std::cout<<"Shift Gradient 0.0"<<endl;
//        }
//        curGrad = newGrad;

//        auto maxForce = maxF();
//        Ms_isminimizationstate = maxForce > GRADTOL;

//#endif
//    }
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
        FFM.computeEnergy(Bead::getDbData().coords.data(), Bead::getDbData().forces.data(), 0.0, true);
#ifdef CUDAACCL
        for(auto strm:CUDAcommon::getCUDAvars().streamvec)
            CUDAcommon::handleerror(cudaStreamSynchronize(*strm));
#endif
        cout << endl;
    }
#if defined(CROSSCHECK) || defined(CUDAACCL)
    cross_checkclass::Aux=false;
#endif
#ifdef CUDAACCL
    cvars = CUDAcommon::getCUDAvars();
    cvars.streamvec.clear();
    CUDAcommon::cudavars = cvars;
#endif

    //final force calculation
    FFM.computeForces(Bead::getDbData().coords.data(), Bead::getDbData().forces.data());
#ifdef ALLSYNC
    cudaDeviceSynchronize();
#endif
#ifdef CUDATIMETRACK
    tbegin = chrono::high_resolution_clock::now();
#endif
#ifdef CUDAACCL
    for(auto strm:CUDAcommon::getCUDAvars().streamvec)
        CUDAcommon::handleerror(cudaStreamSynchronize(*strm));
#endif
#ifdef CUDATIMETRACK
    tend= chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed_run2(tend - tbegin);
    CUDAcommon::cudatime.TveccomputeF.push_back(elapsed_run2.count());
    CUDAcommon::cudatime.TcomputeF += elapsed_run2.count();
#endif

#ifdef CUDAACCL
    FFM.CUDAcopyForces(*Msp1, CUDAcommon::getCUDAvars().gpu_forceAux,CUDAcommon::getCUDAvars().gpu_force);
    //copy back forces and calculate load forces in CPU.
#endif
#ifdef ALLSYNC
    cudaDeviceSynchronize();
#endif
#ifdef SERIAL
    Bead::getDbData().forcesAux = Bead::getDbData().forces;

#endif
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
#ifdef CUDATIMETRACK
    tendII= chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed_runslice4(tendII - tbeginII);
    std::cout<<"Slice time "<<elapsed_runslice4.count()<<endl;
    tbeginII = chrono::high_resolution_clock::now();
#endif
    FFM.computeLoadForces();
    //std::cout<<"End Minimization************"<<endl;
    FFM.cleanupAllForceFields();
#ifdef DETAILEDOUTPUT
    std::cout<<"printing beads & forces"<<endl;
    for(auto b:Bead::getBeads()){
        std::cout<<b->getId()<<" "<< b->coordinate() <<" "
                ""<<b->force() <<endl;
    }
    std::cout<<"printed beads & forces"<<endl;
#endif

#ifdef CUDATIMETRACK
    tendII= chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed_runslice3(tendII - tbeginII);
    std::cout<<"Slice time "<<elapsed_runslice3.count()<<endl;
    tendTot= chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed_runtot(tendTot - tbeginTot);
    std::cout<<"Total Minimization time "<<elapsed_runtot.count()<<endl;
#endif
}
