
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

#include "CGMethod.h"

#include "ForceFieldManager.h"

#include "CGMethodCUDA.h"
#include "Bead.h"
#ifdef CUDAACCL
#ifdef __CUDACC__
#define CUDA_HOSTDEV __host__ __device__
#else
#define CUDA_HOSTDEV
#endif
#ifndef THREADSPERBLOCK
#define THREADSPERBLOCK 512
#endif
#define ARRAY_SIZE 128
//
#include <cuda.h>
#include <cuda_runtime.h>
#include "CUDAcommon.h"
#include "nvToolsExt.h"
#endif
#include <vector>
#include <cmath>
#include <ctime>
#include "Bead.h"
#include <ctime>
#include <cstdlib>
#include "cross_check.h"
#include "nvToolsExt.h"
//
long CGMethod::N = 0;
#ifdef CUDAACCL
void CGMethod::CUDAresetlambda(cudaStream_t stream) {
    resetlambdaCUDA<<<1,1,0, stream>>>(CUDAcommon::getCUDAvars().gpu_lambda);
            CUDAcommon::handleerror(cudaGetLastError(), "resetlambdaCUDA", "CGMethod.cu");
}
void CGMethod::CUDAinitializeLambda(cudaStream_t stream, bool *check_in, bool *check_out, bool *Polaksafestate, int
                                    *gpu_state){
    nvtxRangePushA("CMAXFevstrm");
//    cudaStream_t  s;
//    CUDAcommon::handleerror(cudaStreamCreate(&s));
//    cudaEvent_t  e;
//    CUDAcommon::handleerror(cudaEventCreate(&e));
    nvtxRangePop();
//    nvtxRangePushA("cMAXFv1");
//    maxFCUDA<<<1,1, 0, s>>>(CUDAcommon::getCUDAvars().gpu_forceAux, gpu_nint, gpu_fmax);
//    cudaStreamSynchronize(s);
//    nvtxRangePop();

//    nvtxRangePushA("cMAXFv2");
//    maxFCUDAred<<<1,3, 3*sizeof(double), s>>>(CUDAcommon::getCUDAvars().gpu_forceAux, gpu_nint, gpu_fmax);
//    cudaStreamSynchronize(s);
//    nvtxRangePop();

//    CUDAcommon::handleerror(cudaDeviceSynchronize());
//    std::cout<<"======"<<endl;
//    CUDAcommon::handleerror(cudaEventRecord(e,s));
//    CUDAcommon::handleerror(cudaGetLastError(), "maxFCUDA", "CGMethod.cu");

//    nvtxRangePushA("cMAXFwait");
////    CUDAcommon::handleerror(cudaStreamWaitEvent(stream,e,0));
////    nvtxRangePop();
////    nvtxRangePushA("cMAXFend");
////    CUDAcommon::handleerror(cudaEventDestroy(e));
////    CUDAallFDotFA(s);//
//    CUDAcommon::handleerror(cudaStreamDestroy(s));
//    nvtxRangePop();

    nvtxRangePushA("initialize_lambda");
//    auto gpu_lambda = CUDAcommon::getCUDAvars().gpu_lambda;
    auto gpu_energy = CUDAcommon::getCUDAvars().gpu_energy;
    initializeLambdaCUDA<<<1,1,0, stream>>>(check_in, check_out, g_currentenergy, gpu_energy, gpu_initlambdalocal, gpu_fmax,
            gpu_params, Polaksafestate, gpu_state);
    nvtxRangePop();
//    nvtxRangePushA("init_lambda_sync");
//    CUDAcommon::handleerror(cudaStreamSynchronize (stream));
//    nvtxRangePop();
    CUDAcommon::handleerror(cudaGetLastError(), "initializeLambdaCUDA", "CGMethod.cu");
}

//void CGMethod::getmaxFCUDA(double *gpu_forceAux, int *gpu_nint, double *gpu_fmax) {
//    maxFCUDA<<<1,1>>>(CUDAcommon::getCUDAvars().gpu_forceAux, gpu_nint, gpu_fmax);
//    CUDAcommon::handleerror(cudaGetLastError(), "getmaxFCUDA", "CGMethod.cu");
//}
void CGMethod::CUDAfindLambda(cudaStream_t  stream1, cudaStream_t stream2, cudaEvent_t  event, bool *checkin, bool
        *checkout, bool *gpu_safestate, int *gpu_state) {

    auto gpu_energy = CUDAcommon::getCUDAvars().gpu_energy;
    auto gpu_lambda = CUDAcommon::getCUDAvars().gpu_lambda;
    //update FMAX
//    maxFCUDA<<<1,1, 0, stream>>>(CUDAcommon::getCUDAvars().gpu_forceAux, gpu_nint, gpu_fmax);
////    CUDAcommon::handleerror(cudaDeviceSynchronize());
////    std::cout<<"======"<<endl;
//    CUDAcommon::handleerror(cudaEventRecord(event, stream));
//    CUDAcommon::handleerror(cudaGetLastError(), "maxFCUDA", "CGMethod.cu");
    findLambdaCUDA << < 1, 1, 0, stream1 >> > (gpu_energy, g_currentenergy, gpu_FDotFA, gpu_fmax, gpu_lambda,
            gpu_params, checkin, checkout, gpu_safestate, gpu_state);
    CUDAcommon::handleerror(cudaEventRecord(event, stream1));
    findLambdaCUDA2 << < 1, 1, 0, stream2 >> > (gpu_fmax, gpu_lambda, gpu_params, checkin, checkout, gpu_safestate,
            gpu_state);
    CUDAcommon::handleerror(cudaGetLastError(), "findLambdaCUDA", "CGMethod.cu");
}
//void CGMethod::CUDAprepforbacktracking(cudaStream_t stream, bool *check_in, bool *check_out){
//    nvtxRangePushA("CMAXFevstrm");
//    cudaStream_t  s;
//    CUDAcommon::handleerror(cudaStreamCreate(&s));
//    cudaEvent_t  e;
//    CUDAcommon::handleerror(cudaEventCreate(&e));
//    nvtxRangePop();
//    nvtxRangePushA("cMAXF");
//    maxFCUDA<<<1,1, 0, s>>>(CUDAcommon::getCUDAvars().gpu_forceAux, gpu_nint, gpu_fmax);
//    CUDAcommon::handleerror(cudaEventRecord(e,s));
//    CUDAcommon::handleerror(cudaGetLastError());
//    nvtxRangePop();
//    nvtxRangePushA("cMAXFwait");
//    CUDAcommon::handleerror(cudaStreamWaitEvent(stream,e,0));
//    nvtxRangePop();
//    nvtxRangePushA("cMAXFend");
//    CUDAcommon::handleerror(cudaEventDestroy(e));
//    CUDAcommon::handleerror(cudaStreamDestroy(s));
//    nvtxRangePop();
////    CUDAcommon::handleerror(cudaStreamSynchronize (stream));
//    auto gpu_lambda = CUDAcommon::getCUDAvars().gpu_lambda;
//    auto gpu_energy = CUDAcommon::getCUDAvars().gpu_energy;
//    prepbacktracking<<<1,1,0, stream>>>(check_in, check_out, g_currentenergy, gpu_energy, gpu_lambda, gpu_fmax,
//            gpu_params);
//    CUDAcommon::handleerror(cudaStreamSynchronize (stream));
//    CUDAcommon::handleerror(cudaGetLastError());
//}
//void CGMethod::CUDAprepforsafebacktracking(cudaStream_t stream, bool *check_in, bool *check_out){
//    auto gpu_lambda = CUDAcommon::getCUDAvars().gpu_lambda;
//    auto gpu_energy = CUDAcommon::getCUDAvars().gpu_energy;
//    prepsafebacktracking<<<1,1,0,stream>>>(check_in, check_out, g_currentenergy, gpu_energy, gpu_lambda, gpu_params);
//    CUDAcommon::handleerror(cudaStreamSynchronize (stream));
//    CUDAcommon::handleerror(cudaGetLastError());
//}

void CGMethod::CUDAallFDotF(cudaStream_t stream){
    nvtxRangePushA("COTHR");
    allFADotFCUDA<<<blocksnthreads[0], blocksnthreads[1],0,stream>>>(CUDAcommon::getCUDAvars().gpu_force,
            CUDAcommon::getCUDAvars().gpu_force ,gpu_g, gpu_nint);
    CUDAcommon::handleerror(cudaGetLastError(), "allFADotFCUDA", "CGMethod.cu");
//    addvector<<<1,1,0,stream>>>(gpu_g, gpu_nint, gpu_FDotF);
//    cudaStreamSynchronize(stream);
    addvectorred<<<1,10,10 * sizeof(double),stream>>>(gpu_g, gpu_nint, gpu_FDotF);
//    cudaStreamSynchronize(stream);
    CUDAcommon::handleerror(cudaGetLastError(), "allFADotFCUDA", "CGMethod.cu");
    nvtxRangePop();
}
void CGMethod::CUDAallFADotFA(cudaStream_t stream){
    nvtxRangePushA("COTHR");
    allFADotFCUDA<<<blocksnthreads[0], blocksnthreads[1],0,stream>>>(CUDAcommon::getCUDAvars().gpu_forceAux,
            CUDAcommon::getCUDAvars().gpu_forceAux ,gpu_g, gpu_nint);
    CUDAcommon::handleerror(cudaGetLastError(), "allFADotFCUDA", "CGMethod.cu");
//    addvector<<<1,1,0,stream>>>(gpu_g, gpu_nint, gpu_FADotFA);
//    cudaStreamSynchronize(stream);
    addvectorred<<<1,10,10 * sizeof(double),stream>>>(gpu_g, gpu_nint, gpu_FADotFA);
//    cudaStreamSynchronize(stream);
    CUDAcommon::handleerror(cudaGetLastError(), "allFADotFCUDA", "CGMethod.cu");
    nvtxRangePop();
}
void CGMethod::CUDAallFADotFAP(cudaStream_t stream){
    nvtxRangePushA("COTHR");
    allFADotFCUDA<<<blocksnthreads[0], blocksnthreads[1],0,stream>>>(CUDAcommon::getCUDAvars().gpu_forceAux,
            CUDAcommon::getCUDAvars().gpu_forceAuxP ,gpu_g, gpu_nint);
    CUDAcommon::handleerror(cudaGetLastError(), "allFADotFCUDA", "CGMethod.cu");
//    addvector<<<1,1,0,stream>>>(gpu_g, gpu_nint, gpu_FADotFAP);
//    cudaStreamSynchronize(stream);
    addvectorred<<<1,10,10 * sizeof(double),stream>>>(gpu_g, gpu_nint, gpu_FADotFAP);
//    cudaStreamSynchronize(stream);
    CUDAcommon::handleerror(cudaGetLastError(), "allFADotFCUDA", "CGMethod.cu");
    nvtxRangePop();
}
void CGMethod::CUDAallFDotFA(cudaStream_t stream){
    nvtxRangePushA("COTHR");
    allFADotFCUDA<<<blocksnthreads[0], blocksnthreads[1],0,stream>>>(CUDAcommon::getCUDAvars().gpu_force,
            CUDAcommon::getCUDAvars().gpu_forceAux ,gpu_g, gpu_nint);
    CUDAcommon::handleerror(cudaGetLastError(), "allFADotFCUDA", "CGMethod.cu");
//    addvector<<<1,1,0,stream>>>(gpu_g, gpu_nint, gpu_FDotFA);
//    cudaStreamSynchronize(stream);
    addvectorred<<<1,10,10 * sizeof(double),stream>>>(gpu_g, gpu_nint, gpu_FDotFA);
//    cudaStreamSynchronize(stream);
    CUDAcommon::handleerror(cudaGetLastError(), "allFADotFCUDA", "CGMethod.cu");
    nvtxRangePop();
}

void CGMethod::CUDAshiftGradient(cudaStream_t stream, bool *Mcheckin) {
    shiftGradientCUDA<<<blocksnthreads[0], blocksnthreads[1],0, stream>>>(CUDAcommon::getCUDAvars().gpu_force,
            CUDAcommon::getCUDAvars().gpu_forceAux, gpu_nint, gpu_FADotFA, gpu_FADotFAP, gpu_FDotF, Mcheckin);
}

void CGMethod::CUDAshiftGradientifSafe(cudaStream_t stream, bool *Mcheckin, bool *Scheckin){
    shiftGradientCUDAifsafe<<<blocksnthreads[0], blocksnthreads[1],0, stream>>>(CUDAcommon::getCUDAvars().gpu_force, CUDAcommon::getCUDAvars().gpu_forceAux, gpu_nint,
                            Mcheckin, Scheckin);
}

//void CGMethod::CUDAgetPolakvars(bool calc_safestate,cudaStream_t streamcalc, double* gpu_GRADTOL, bool *gminstatein,
//                                    bool *gminstateout, bool *gsafestateout, volatile bool *cminstate){
////    state[0] = false;
////    state[1] = false;
//    if(cminstate[0] == true) {
//        nvtxRangePushA("cMAXF");
////        maxFCUDA << < 1, 1, 0, streamcalc >> > (CUDAcommon::getCUDAvars().gpu_forceAux, gpu_nint, gpu_fmax);
//        maxFCUDAred<<<1,3, 3*sizeof(double), streamcalc>>>(CUDAcommon::getCUDAvars().gpu_forceAux, gpu_nint, gpu_fmax);
////        CUDAcommon::handleerror(cudaDeviceSynchronize());
////        std::cout<<"======"<<endl;
//        CUDAcommon::handleerror(cudaGetLastError(), "maxFCUDA", "CGMethod.cu");
//        nvtxRangePop();
//        getminimizestateCUDA << < 1, 1, 0, streamcalc >> > (gpu_fmax, gpu_GRADTOL, gminstatein, gminstateout);
//        CUDAcommon::handleerror(cudaGetLastError(), "getminimizestateCUDA", "CGMethod.cu");
////        CUDAcommon::handleerror(cudaStreamSynchronize(streamcalc));
//    }
//    if(calc_safestate){
//        CUDAallFDotFA(streamcalc);
//        getsafestateCUDA<<<1,1,0,streamcalc>>>(gpu_FDotFA, gpu_FDotF, gpu_FADotFA, gsafestateout);
//        CUDAcommon::handleerror(cudaGetLastError(), "getsafestateCUDA", "CGMethod.cu");
//    }
//}

void CGMethod::CUDAgetPolakvars(cudaStream_t streamcalc, double* gpu_GRADTOL, bool *gminstatein,
                                bool *gminstateout, volatile bool *cminstate){
//    state[0] = false;
//    state[1] = false;
    if(cminstate[0] == true) {
        nvtxRangePushA("cMAXF");
//        maxFCUDA << < 1, 1, 0, streamcalc >> > (CUDAcommon::getCUDAvars().gpu_forceAux, gpu_nint, gpu_fmax);
//        maxFCUDAred<<<1,3, 3*sizeof(double), streamcalc>>>(CUDAcommon::getCUDAvars().gpu_forceAux, gpu_nint, gpu_fmax);
//        cudaStreamSynchronize(streamcalc);
        maxFCUDAredv2<<<1,10, 10*sizeof(double), streamcalc>>>(CUDAcommon::getCUDAvars().gpu_forceAux, gpu_nint,
                gpu_fmax);
//        cudaStreamSynchronize(streamcalc);
//        CUDAcommon::handleerror(cudaDeviceSynchronize());
//        std::cout<<"======"<<endl;
        CUDAcommon::handleerror(cudaGetLastError(), "maxFCUDA", "CGMethod.cu");
        nvtxRangePop();
        getminimizestateCUDA << < 1, 1, 0, streamcalc >> > (gpu_fmax, gpu_GRADTOL, gminstatein, gminstateout);
        CUDAcommon::handleerror(cudaGetLastError(), "getminimizestateCUDA", "CGMethod.cu");
//        CUDAcommon::handleerror(cudaStreamSynchronize(streamcalc));
    }
//    if(calc_safestate){
//        CUDAallFDotFA(streamcalc);
//        getsafestateCUDA<<<1,1,0,streamcalc>>>(gpu_FDotFA, gpu_FDotF, gpu_FADotFA, gsafestateout);
//        CUDAcommon::handleerror(cudaGetLastError(), "getsafestateCUDA", "CGMethod.cu");
//    }
}

void CGMethod::CUDAgetPolakvars2(cudaStream_t streamcalc, bool *gsafestateout){
        CUDAallFDotFA(streamcalc);
        getsafestateCUDA<<<1,1,0,streamcalc>>>(gpu_FDotFA, gpu_FDotF, gpu_FADotFA, gsafestateout);
        CUDAcommon::handleerror(cudaGetLastError(), "getsafestateCUDA", "CGMethod.cu");

}

void CGMethod::CUDAmoveBeads(cudaStream_t stream, bool *gpu_checkin){
    double *gpu_lambda = CUDAcommon::getCUDAvars().gpu_lambda;
    double *gpu_coord = CUDAcommon::getCUDAvars().gpu_coord;
    double *gpu_force = CUDAcommon::getCUDAvars().gpu_force;

    nvtxRangePushA("CMVB");

    moveBeadsCUDA<<<blocksnthreads[0], blocksnthreads[1],0, stream>>>(gpu_coord, gpu_force, gpu_lambda, gpu_nint,
            gpu_checkin);
    nvtxRangePop();
    CUDAcommon::handleerror(cudaGetLastError(),"moveBeadsCUDA", "CGMethod.cu");
}

void CGMethod::CUDAinitializePolak(cudaStream_t stream, bool *minstatein, bool *minstateout, bool *safestatein, bool
        *safestateout){
    initializePolak<<<1,1,0,stream>>>(minstatein, minstateout, safestatein, safestateout);
}

double CGMethod::gpuFDotF(double *f1,double *f2){

    allFADotFCUDA<<<blocksnthreads[0], blocksnthreads[1]>>>(f1, f2 ,gpu_g, gpu_nint);
    CUDAcommon::handleerror(cudaGetLastError(),"allFADotFCUDA", "CGMethod.cu");
//    addvector<<<1,1>>>(gpu_g, gpu_nint, gSum);
    addvectorred<<<1,10,10 * sizeof(double)>>>(gpu_g, gpu_nint, gSum);
    CUDAcommon::handleerror(cudaGetLastError(),"allFADotFCUDA", "CGMethod.cu");

//    CUDAcommon::handleerror( cudaPeekAtLastError() );
//    CUDAcommon::handleerror(cudaDeviceSynchronize());

    double g[1];
    CUDAcommon::handleerror(cudaMemcpy(g, gSum, sizeof(double),
                                       cudaMemcpyDeviceToHost));


//    double g[N/3];
//    CUDAcommon::handleerror(cudaMemcpy(g, gpu_g, N/3 * sizeof(double),
//                                       cudaMemcpyDeviceToHost));
//    CUDAcommon::handleerror(cudaFree(gpu_g));
//    double sum=0.0;
//    for(auto i=0;i<N/3;i++)
//        sum+=g[i];
    return g[0];
}
#endif
double CGMethod::allFDotF()
{
    nvtxRangePushA("SOTHR");
    double g = 0;
    for(int i = 0; i < N; i++)
        g += force[i] * force[i];
    nvtxRangePop();
    return g;
}

double CGMethod::allFADotFA()
{
    nvtxRangePushA("SOTHR");
    double g = 0;
    for(int i = 0; i < N; i++)
        g += forceAux[i] * forceAux[i];
    nvtxRangePop();
//#ifdef CUDAACCL
//    nvtxRangePushA("COTHR");
//    auto g_cuda = gpuFDotF(CUDAcommon::getCUDAvars().gpu_forceAux,CUDAcommon::getCUDAvars().gpu_forceAux);
//    nvtxRangePop();
//#endif
//    if(g>1000000.0){
//        if(abs(g-g_cuda)/abs(g) > 0.001){
//            std::cout << g << " " << g_cuda << endl;
//            std::cout << "Precison mismatch FADotFA " << abs(g - g_cuda) << endl;
//        }
//
//    }
//    else if(abs(g-g_cuda)>1/100000000.0) {
//        std::cout << g << " " << g_cuda << endl;
//        std::cout << "Precison mismatch FADotFA " << abs(g - g_cuda) << endl;
//    }
    return g;
}

double CGMethod::allFADotFAP()
{
    nvtxRangePushA("SOTHR");
    double g = 0;
    for(int i = 0; i < N; i++)
        g += forceAux[i] * forceAuxPrev[i];
    nvtxRangePop();
    return g;
}

double CGMethod::allFDotFA()
{
    nvtxRangePushA("SOTHR");
    double g = 0;
    for(int i = 0; i < N; i++)
        g += force[i] * forceAux[i];
    nvtxRangePop();
//#ifdef CUDAACCL
//    nvtxRangePushA("COTHR");
//    auto g_cuda = gpuFDotF(CUDAcommon::getCUDAvars().gpu_force,CUDAcommon::getCUDAvars().gpu_forceAux);
//    nvtxRangePop();
//#endif
//    if(g>1000000.0){
//        if(abs(g-g_cuda)/abs(g) > 0.001){
//            std::cout << g << " " << g_cuda << endl;
//            std::cout << "Precison mismatch FDotFA " << abs(g - g_cuda) << endl;
//        }
//
//    }
//    else if(abs(g-g_cuda)>1/100000000.0) {
//        std::cout << "Precison mismatch FDotFA " << abs(g - g_cuda) << endl;
//        std::cout << g << " " << g_cuda << endl;
//
//    }
    return g;
}

double CGMethod::maxF() {
    nvtxRangePushA("SMAXF");
    double maxF = 0;
    double mag;

    for(int i = 0; i < N; i++) {
        mag = sqrt(forceAux[i]*forceAux[i]);
        if(mag > maxF) maxF = mag;
    }
    nvtxRangePop();
    return maxF;
}

Bead* CGMethod::maxBead() {

    double maxF = 0.0;
    double currentF;
    long index = 0;
#ifndef CUDAACCL
    for (int i = 0; i < N; i++) {

        currentF = forceAux[i] * forceAux[i];
        if(currentF > maxF) {
            index = i;
            maxF = currentF;
        }
    }
#endif
#ifdef CUDAACCL
    double F_i[N];
    double gmaxF = 0.0;
    CUDAcommon::handleerror(cudaDeviceSynchronize());
    CUDAcommon::handleerror(cudaMemcpy(F_i, CUDAcommon::getCUDAvars().gpu_forceAux, N *
                                                                                 sizeof(double), cudaMemcpyDeviceToHost));
    double gcurrentF;
//    long gindex = 0;

    for (int i = 0; i < N; i++) {

        gcurrentF = F_i[i] * F_i[i];
        if(gcurrentF > gmaxF) {
            index = i;
            gmaxF = gcurrentF;
//            std::cout<<gcurrentF<<" "<<forceAux[i] * forceAux[i]<<endl;
        }
    }
//    if(gindex!=index)
//        std::cout<<N<<endl;
//        std::cout<<"CPU and GPU codes do not point to same bead with maxF."<<endl;
#endif
    return Bead::getBeads()[index];
}

void CGMethod::moveBeads(double d)
{
    ///<NOTE: Ignores static beads for now.
    //if(!b->getstaticstate())
    nvtxRangePushA("SMVB");
    for (int i = 0; i < N; i++)
        coord[i] = coord[i] + d * force[i];
    nvtxRangePop();
}

void CGMethod::shiftGradient(double d)
{
    for (int i = 0; i < N; i ++)
        force[i] = forceAux[i] + d * force[i];
}

void CGMethod::printForces()
{
    cout << "Print Forces" << endl;
    for(auto b: Bead::getBeads()) {

        for (int i = 0; i<3; i++)
            cout << b->coordinate[i] << "  "<<
                 b->force[i] <<"  "<<b->forceAux[i]<<endl;
    }
    cout << "End of Print Forces" << endl;
}

void CGMethod::startMinimization() {


    //COPY BEAD DATA
    N = 3 * Bead::getBeads().size();
//    std::cout<<3 * Bead::getBeads().size()<<endl;
    allocate(N);

    //coord management
    long i = 0;
    long index = 0;
    for(auto b: Bead::getBeads()) {

        //set bead index
        b->_dbIndex = i;

        //flatten indices
        index = 3 * i;
        coord[index] = b->coordinate[0];
        coord[index + 1] = b->coordinate[1];
        coord[index + 2] = b->coordinate[2];

        b->coordinateP = b->coordinate;
        i++;
    }

#ifdef CUDAACCL

    int nDevices;
//    cudaDeviceProp prop;
    cudaGetDeviceCount(&nDevices);
    if(nDevices>1){
        cout<<"Code not configured for multiple devices. Exiting..."<<endl;
        exit(EXIT_FAILURE);
    }
//    for (int i = 0; i < nDevices;  i++) {
//        cudaGetDeviceProperties(&prop, i);
//        printf("Device Number: %d\n", i);
//        printf("  Device name: %s\n", prop.name);
//        printf("  Memory Clock Rate (KHz): %d\n",
//               prop.memoryClockRate);
//        printf("  Memory Bus Width (bits): %d\n",
//               prop.memoryBusWidth);
//        printf("  Peak Memory Bandwidth (GB/s): %f\n\n",
//               2.0*prop.memoryClockRate*(prop.memoryBusWidth/8)/1.0e6);
//        printf("  Total amount of shared memory per block:       %u bytes\n", prop.sharedMemPerBlock);
//            printf("  Total amount of global memory:                 %llu bytes\n", (unsigned long long) prop
//                    .totalGlobalMem);
//    }

    nvtxRangePushA("CSM");

    double f[N];
    for(auto iter=0;i<N;i++)
        f[iter]=0.0;
    double* gpu_coord;
    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_coord, N*sizeof(double)));
    double* gpu_lambda;
    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_lambda, sizeof(double)));
    double* gpu_force;
    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_force, N*sizeof(double)));
    double* gpu_forceAux;
    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_forceAux, N*sizeof(double)));
    double* gpu_forceAuxP;
    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_forceAuxP, N*sizeof(double)));
    double* gpu_energy;
    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_energy, sizeof(double)));
    bool* gpu_btstate;
    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_btstate, sizeof(bool)));

    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_initlambdalocal, sizeof(double)));

    CUDAcommon::handleerror(cudaMalloc((void **)&gpu_fmax, sizeof(double)));
    CUDAcommon::handleerror(cudaMalloc((void **)&g_currentenergy, sizeof(double)));
    CUDAcommon::handleerror(cudaMalloc((void **)&gpu_FDotF, sizeof(double)));
    CUDAcommon::handleerror(cudaMalloc((void **)&gpu_FADotFA, sizeof(double)));
    CUDAcommon::handleerror(cudaMalloc((void **)&gpu_FADotFAP, sizeof(double)));
    CUDAcommon::handleerror(cudaMalloc((void **)&gpu_FDotFA, sizeof(double)));

    CUDAcommon::handleerror(cudaHostAlloc((void**)&convergencecheck, 3 * sizeof(bool), cudaHostAllocMapped));
    CUDAcommon::handleerror(cudaHostGetDevicePointer(&gpu_convergencecheck, convergencecheck, 0));

    //PING PONG
    cudaMalloc(&g_stop1, sizeof(bool));
    cudaMalloc(&g_stop2, sizeof(bool));
    cudaHostAlloc(&h_stop, sizeof(bool), cudaHostAllocDefault);
    //@

//    CUDAcommon::handleerror(cudaHostAlloc((void**)&convergencecheck, 3 * sizeof(bool), cudaHostAllocDefault));
//    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_convergencecheck, 3 * sizeof(bool)));

//    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_lambda, sizeof(double))); REPEAT.
    CUDAcommon::handleerror(cudaMemcpy(gpu_coord, coord, N*sizeof(double), cudaMemcpyHostToDevice));
    CUDAcommon::handleerror(cudaMemcpy(gpu_force, f, N*sizeof(double), cudaMemcpyHostToDevice));
    CUDAcommon::handleerror(cudaMemcpy(gpu_forceAux, f, N*sizeof(double), cudaMemcpyHostToDevice));
    CUDAcommon::handleerror(cudaMemcpy(gpu_forceAuxP, f, N*sizeof(double), cudaMemcpyHostToDevice));
    bool dummy[1];dummy[0] = true;
    CUDAcommon::handleerror(cudaMemcpy(gpu_btstate, dummy, sizeof(bool), cudaMemcpyHostToDevice));
    int *gculpritID;
    char *gculpritFF;
    char *gculpritinteraction;
    int *culpritID;
    char *culpritFF;
    char *culpritinteraction;
    CUDAcommon::handleerror(cudaHostAlloc((void**)&culpritID, 4 * sizeof(int), cudaHostAllocMapped));
    CUDAcommon::handleerror(cudaHostAlloc((void**)&culpritFF, 100*sizeof(char), cudaHostAllocMapped));
    CUDAcommon::handleerror(cudaHostAlloc((void**)&culpritinteraction, 100*sizeof(char), cudaHostAllocMapped));
    CUDAcommon::handleerror(cudaHostGetDevicePointer(&gculpritID, culpritID, 0));
    CUDAcommon::handleerror(cudaHostGetDevicePointer(&gculpritFF, culpritFF, 0));
    CUDAcommon::handleerror(cudaHostGetDevicePointer(&gculpritinteraction, culpritinteraction, 0));
//    CUDAcommon::handleerror(cudaMalloc((void **) &gculpritID, sizeof(int)));
//    CUDAcommon::handleerror(cudaMalloc((void **) &gculpritFF, 11*sizeof(char)));
//    char a[] = "FilamentFF";
//    CUDAcommon::handleerror(cudaMemcpy(gculpritFF, a, 100 * sizeof(char), cudaMemcpyHostToDevice));
//    CUDAcommon::handleerror(cudaMalloc((void **) &gculpritinteraction, 100*sizeof(char)));

    CUDAvars cvars=CUDAcommon::getCUDAvars();
    cvars.gpu_coord=gpu_coord;
    cvars.gpu_lambda=gpu_lambda;
    cvars.gpu_forceAux = gpu_forceAux;
    cvars.gpu_force = gpu_force;
    cvars.gpu_forceAuxP = gpu_forceAuxP;
    cvars.gpu_energy = gpu_energy;
    cvars.gculpritID = gculpritID;
    cvars.culpritID = culpritID;
    cvars.gculpritinteraction = gculpritinteraction;
    cvars.gculpritFF = gculpritFF;
    cvars.culpritinteraction = culpritinteraction;
    cvars.culpritFF = culpritFF;
    cvars.gpu_btstate = gpu_btstate;
    CUDAcommon::cudavars=cvars;
//SET CERTAIN GPU PARAMETERS SET FOR EASY ACCESS DURING MINIMIZATION.
    int nint[1]; nint[0]=CGMethod::N/3;
    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_nint, sizeof(int)));
    CUDAcommon::handleerror(cudaMemcpy(gpu_nint, nint, sizeof(int), cudaMemcpyHostToDevice));
    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_state, sizeof(int)));
    blocksnthreads.push_back(CGMethod::N/(3*THREADSPERBLOCK) + 1);
    if(blocksnthreads[0]==1) blocksnthreads.push_back(CGMethod::N/3);
    else blocksnthreads.push_back(THREADSPERBLOCK);
    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_g, N/3 * sizeof(double)));
    CUDAcommon::handleerror(cudaMalloc((void **) &gSum, sizeof(double)));
    nvtxRangePop();

//    cvars.gpu_globalMem = prop.totalGlobalMem;
//    cvars.gpu_sharedMem = prop.sharedMemPerBlock;
//    double a;
//    std::cout<<cvars.gpu_globalMem<<" "<<cvars.gpu_sharedMem<<" "<<sizeof(a)<<endl;
//
//    double ccoord[N];
//    cudaMemcpy(ccoord, gpu_coord, N*sizeof(double), cudaMemcpyDeviceToHost);
//    for(auto i=0;i<N;i++)
//        std::cout<<ccoord[i]<<" "<<coord[i]<<endl;

//    vector<double> c2;c2.push_back(273.14);c2.push_back(273.14);
//    double c2[2];
//    c2[0]=10.234;c2[1]=20.234;
//    double *gpu_coord2;
//    cudaMalloc((void **) &gpu_coord2, 2*sizeof(double));
//    cudaMemcpy(gpu_coord2, c2, 2*sizeof(double), cudaMemcpyHostToDevice);
//
//    double cc[2];
//    cudaMemcpy(cc, gpu_coord2, 2*sizeof(double), cudaMemcpyDeviceToHost);
//    std::cout<<cc[0]<<" "<<cc[1]<<endl;
//    cudaFree(gpu_coord2);
//    cudaFree(gpu_coord);
#endif
}

void CGMethod::endMinimization() {
#ifdef CUDAACCL
    nvtxRangePushA("CEM");
    CUDAcommon::handleerror(cudaMemcpy(coord, CUDAcommon::getCUDAvars().gpu_coord, N *
                            sizeof(double), cudaMemcpyDeviceToHost));
    nvtxRangePop();
    #endif
    ///RECOPY BEAD DATA
    //coord management
    long i = 0;
    long index = 0;
    for(auto b: Bead::getBeads()) {

        //flatten indices
        index = 3 * i;
        b->coordinate[0] = coord[index];
        b->coordinate[1] = coord[index + 1];
        b->coordinate[2] = coord[index + 2];
        i++;
    }

    deallocate();
#ifdef CUDAACCL
    nvtxRangePushA("CEM");
    CUDAcommon::handleerror(cudaFree(CUDAcommon::getCUDAvars().gpu_coord));
    CUDAcommon::handleerror(cudaFree(CUDAcommon::getCUDAvars().gpu_force));
    CUDAcommon::handleerror(cudaFree(CUDAcommon::getCUDAvars().gpu_forceAux));
    CUDAcommon::handleerror(cudaFree(CUDAcommon::getCUDAvars().gpu_forceAuxP));
    CUDAcommon::handleerror(cudaFree(CUDAcommon::getCUDAvars().gpu_lambda));
    CUDAcommon::handleerror(cudaFree(CUDAcommon::getCUDAvars().gpu_energy));
//    CUDAcommon::handleerror(cudaFree(CUDAcommon::getCUDAvars().gculpritID));

    CUDAcommon::handleerror(cudaFree(CUDAcommon::getCUDAvars().gpu_btstate));
    CUDAcommon::handleerror(cudaFree(gpu_initlambdalocal));
//    CUDAcommon::handleerror(cudaFreeHost(CUDAcommon::getCUDAvars().culpritFF));
    CUDAcommon::handleerror(cudaFreeHost(CUDAcommon::getCUDAvars().culpritID));
    CUDAcommon::handleerror(cudaFreeHost(CUDAcommon::getCUDAvars().culpritFF));
    CUDAcommon::handleerror(cudaFreeHost(CUDAcommon::getCUDAvars().culpritinteraction));
    CUDAcommon::handleerror(cudaFree(gpu_g));
    CUDAcommon::handleerror(cudaFree(gSum));
    CUDAcommon::handleerror(cudaFree(gpu_fmax));
    CUDAcommon::handleerror(cudaFree(gpu_FDotF));
    CUDAcommon::handleerror(cudaFree(gpu_FADotFA));
    CUDAcommon::handleerror(cudaFree(gpu_FADotFAP));
    CUDAcommon::handleerror(cudaFree(gpu_FDotFA));
    CUDAcommon::handleerror(cudaFreeHost(convergencecheck));
    CUDAcommon::handleerror(cudaFree(g_currentenergy));
    //PING PONG SAFEBACKTRACKING AND BACKTRACKING
    CUDAcommon::handleerror(cudaFree(g_stop1));
    CUDAcommon::handleerror(cudaFree(g_stop2));
    CUDAcommon::handleerror(cudaFreeHost(h_stop));
    //@
//    CUDAcommon::handleerror(cudaFree(gpu_convergencecheck));


    CUDAcommon::handleerror(cudaFree(gpu_nint));
    CUDAcommon::handleerror(cudaFree(gpu_state));
    blocksnthreads.clear();
    nvtxRangePop();
    //TODO cross check later
//    CUDAcommon::handleerror(cudaFree(CUDAcommon::getCUDAvars().motorparams));

//    CUDAcommon::getCUDAvars().gpu_coord = NULL;
//    CUDAcommon::getCUDAvars().gpu_force = NULL;
//    CUDAcommon::getCUDAvars().gpu_forceAux = NULL;
//    CUDAcommon::getCUDAvars().gpu_lambda = NULL;
#endif
}

double CGMethod::backtrackingLineSearch(ForceFieldManager& FFM, double MAXDIST,
                                        double LAMBDAMAX, bool *gpu_safestate) {
    double lambda;
    sconvergencecheck = true;
#ifndef CUDAACCL //SERIAL
    sconvergencecheck = false;
    cconvergencecheck = new bool[1];
    cconvergencecheck[0] = true;
#endif
#ifdef CUDAACCL
    nvtxRangePushA("Evcreate");
    CUDAcommon::handleerror(cudaStreamCreate(&s1));
    CUDAcommon::handleerror(cudaStreamCreate(&s2));
    CUDAcommon::handleerror(cudaStreamCreate(&s3));
    CUDAcommon::handleerror(cudaEventCreate(&e1));
    CUDAcommon::handleerror(cudaEventCreate(&e2));
    nvtxRangePop();
    sp1 = &s1;
    sp2 = &s2;
    ep1 = &e1;
    ep2 = &e2;
    g_s1 = g_stop1;
    g_s2 = g_stop2;
    //prep for backtracking.
    if(gpu_params ==NULL){
        double params[5];
        params[0] = BACKTRACKSLOPE;
        params[1] = LAMBDAREDUCE;
        params[2] = LAMBDATOL;
        params[3] = LAMBDAMAX;
        params[4] = MAXDIST;
        CUDAcommon::handleerror(cudaMalloc((void **) &gpu_params, 5 * sizeof(double)));
        CUDAcommon::handleerror(cudaMemcpy(gpu_params, params, 5 * sizeof(double),
                                           cudaMemcpyHostToDevice));
    }
//    //initialize lambda search
//    CUDAinitializeLambda(*sp1, g_s1, g_s2, gpu_safestate);
//    CUDAcommon::handleerror(cudaEventRecord(*ep1, *sp1));
//    //check if converged.
//    nvtxRangePushA("backConv");
//    CUDAcommon::handleerror(cudaStreamWaitEvent(s3, *ep1, 0));
//    CUDAcommon::handleerror(cudaMemcpyAsync(h_stop, g_s2, sizeof(bool), cudaMemcpyDeviceToHost, s3));
//    nvtxRangePop();

    nvtxRangePushA("reset_lambda_sync");
    CUDAresetlambda(*sp1);//set lambda to zero.
    cudaEvent_t  e;
    CUDAcommon::handleerror(cudaEventCreate(&e));
    CUDAcommon::handleerror(cudaEventRecord(e, *sp1));
//    cudaEventSynchronize(e);
//    cudaStreamSynchronize(*sp1);
//    cudaEventDestroy(e);
    nvtxRangePop();
    auto cvars = CUDAcommon::getCUDAvars();
    cvars.streamvec.clear();
//    cvars.event = &e;
    CUDAcommon::cudavars = cvars;
    //initialize lambda search
    nvtxRangePushA("cuda_init_lambda");
    CUDAinitializeLambda(*sp1, g_s1, g_s2, gpu_safestate, gpu_state);
//    CUDAcommon::handleerror(cudaEventRecord(*ep1, *sp1));
//    cudaStreamSynchronize(*sp1);
    nvtxRangePop();
#else
    double f = maxF();
    //return zero if no forces
    if(f == 0.0){lambda = 0.0; sconvergencecheck = true;}
    //calculate first lambda
    lambda = min(LAMBDAMAX, MAXDIST / f);
#endif
//    std::cout<<"back 0"<<endl;
    nvtxRangePushA("Energy 0");
    double currentEnergy = FFM.computeEnergy(coord, force, 0.0);
    nvtxRangePop();

#ifdef CUDAACCL
    //wait for energies to be calculated
    cudaStream_t stream_bt;
    CUDAcommon::handleerror(cudaStreamCreate(&stream_bt),"find lambda", "CGMethod.cu");
    nvtxRangePushA("backConvsync");
    for(auto strm:CUDAcommon::getCUDAvars().streamvec)
        CUDAcommon::handleerror(cudaStreamSynchronize(*strm),"backVonvSync","CGMethod.cu");
    nvtxRangePop();
    nvtxRangePushA("set_cur_E");
    cudaStreamSynchronize(*sp1);
    setcurrentenergy<<<1,1,0,*sp1>>>(CUDAcommon::getCUDAvars().gpu_energy, g_currentenergy, CUDAcommon::getCUDAvars()
            .gpu_lambda, gpu_initlambdalocal);
    CUDAcommon::handleerror(cudaGetLastError(),"setcurrentenergy", "CGMethod.cu");
    cudaStreamSynchronize(*sp1);
    nvtxRangePop();

    //check if converged.
    nvtxRangePushA("backConv");
    CUDAcommon::handleerror(cudaStreamWaitEvent(s3, *ep1, 0));
//    CUDAcommon::handleerror(cudaEventRecord(*CUDAcommon::getCUDAvars().event, *sp1));
    CUDAcommon::handleerror(cudaMemcpyAsync(h_stop, g_s2, sizeof(bool), cudaMemcpyDeviceToHost, s3));
    nvtxRangePop();
//    CUDAcommon::handleerror(cudaStreamSynchronize (*sp1)); CHECK IF NEEDED
    h_stop[0] = false;
    cconvergencecheck = h_stop;
#endif
int iter = 0;
    while(!(cconvergencecheck[0])||!(sconvergencecheck)) {

        iter++;

#ifdef CUDAACCL
        nvtxRangePushA("While wait");
        CUDAcommon::handleerror(cudaStreamWaitEvent(*sp2, *ep1, 0));
        CUDAcommon::handleerror(cudaStreamSynchronize(*sp2));
        nvtxRangePop();
        //ping pong swap
        sps = sp1;
        sp1 = sp2;
        sp2 = sps;
        eps = ep1;
        ep1 = ep2;
        ep2 = eps;
//        g_bs = g_b1;
//        g_b1 = g_b2;
//        g_b2 = g_bs;
        g_ss = g_s1;
        g_s1 = g_s2;
        g_s2 = g_ss;

        auto cvars = CUDAcommon::getCUDAvars();
        cvars.streamvec.clear();
//        cvars.event = ep1;
        CUDAcommon::cudavars = cvars;
#endif
//        std::cout<<"back Z"<<endl;
        nvtxRangePushA("Energy Z");
        //let each forcefield calculate energy IFF conv state = false. That will help them avoid unnecessary iterations.
        //let each forcefield also add energies to two different energy variables.
        double energyLambda = FFM.computeEnergy(coord, force, lambda);
        nvtxRangePop();

#ifdef CUDAACCL
        //wait for energies to be calculated
//        std::cout<<"Total energy streams "<<CUDAcommon::getCUDAvars().streamvec.size()<<endl;
        nvtxRangePushA("backConvsync");
        for(auto strm:CUDAcommon::getCUDAvars().streamvec)
            CUDAcommon::handleerror(cudaStreamSynchronize(*strm),"backConvsync","CGMethod.cu");
        nvtxRangePop();
//        CUDAcommon::handleerror(cudaDeviceSynchronize());
        if(!(cconvergencecheck[0])){
            std::cout<<"----------"<<endl;
            nvtxRangePushA("While wait2");
            CUDAcommon::handleerror(cudaStreamSynchronize(stream_bt));
            nvtxRangePop();
            nvtxRangePushA("CUDAfindlambda");
            CUDAfindLambda(*sp1, stream_bt, *ep1, g_s1, g_s2, gpu_safestate, gpu_state);
                    CUDAcommon::handleerror(cudaStreamSynchronize(*sp1));
            CUDAcommon::handleerror(cudaStreamSynchronize(stream_bt));
            auto cvars = CUDAcommon::getCUDAvars();
//            cvars.event = ep1;
            CUDAcommon::cudavars = cvars;
            nvtxRangePop();
//            nvtxRangePushA("CUDAfindlambda_sync");
//            CUDAcommon::handleerror(cudaEventSynchronize(*ep1));
//            nvtxRangePop();
            //check if converged.
            nvtxRangePushA("backConv");
            if(cconvergencecheck[0]  == false){
                CUDAcommon::handleerror(cudaStreamWaitEvent(s3, *ep1, 0));
                CUDAcommon::handleerror(cudaMemcpyAsync(h_stop, g_s2, sizeof(bool), cudaMemcpyDeviceToHost, s3));
            }
            nvtxRangePop();
        }
#else

        if(!(sconvergencecheck)){
                double idealEnergyChange = -BACKTRACKSLOPE * lambda * allFDotFA();
                double energyChange = energyLambda - currentEnergy;

                //return if ok
                if(energyChange <= idealEnergyChange) {
//            double cudalambda[1];
//            CUDAcommon::handleerror(cudaMemcpy(cudalambda, CUDAcommon::getCUDAvars().gpu_lambda, sizeof(double),
//                                               cudaMemcpyDeviceToHost));
//            std::cout<<lambda<<" "<<cudalambda[0]<<" "<<convergencecheck[0]<< endl;
                    sconvergencecheck = true;}
                else
                //reduce lambda
                lambda *= LAMBDAREDUCE;

                if(lambda <= 0.0 || lambda <= LAMBDATOL) {
                    sconvergencecheck = true;
                    lambda = 0.0;
//            nvtxRangePushA("CLCP");
//            CUDAcommon::handleerror(cudaMemcpy(CUDAcommon::getCUDAvars().gpu_lambda, &lambda, sizeof(double),
//                                               cudaMemcpyHostToDevice));
//            nvtxRangePop();
//            double cudalambda[1];
//            CUDAcommon::handleerror(cudaMemcpy(cudalambda, CUDAcommon::getCUDAvars().gpu_lambda, sizeof(double),
//                                               cudaMemcpyDeviceToHost));
//            std::cout<<lambda<<" "<<cudalambda[0]<<" "<<convergencecheck[0]<< endl;
                }
        }
#endif
    }
#ifdef CUDAACCL
    nvtxRangePushA("Evsync");
    correctlambdaCUDA<<<1,1,0, stream_bt>>>(CUDAcommon::getCUDAvars().gpu_lambda, gpu_state, gpu_params);
    CUDAcommon::handleerror(cudaStreamSynchronize(stream_bt));
    CUDAcommon::handleerror(cudaStreamSynchronize(s1));
    CUDAcommon::handleerror(cudaStreamSynchronize(s2));
    CUDAcommon::handleerror(cudaStreamSynchronize(s3));
    nvtxRangePop();

    nvtxRangePushA("EvDESTROY");
    CUDAcommon::handleerror(cudaStreamDestroy(s1));
    CUDAcommon::handleerror(cudaStreamDestroy(s2));
    CUDAcommon::handleerror(cudaStreamDestroy(s3));
    CUDAcommon::handleerror(cudaEventDestroy(e1));
    CUDAcommon::handleerror(cudaEventDestroy(e2));
    nvtxRangePop();
#else
    delete cconvergencecheck;
#endif
    std::cout<<"lambda determined in "<<iter<<endl;
//synchronize streams
    if(cconvergencecheck[0]||sconvergencecheck)
        return lambda;

}

double CGMethod::safeBacktrackingLineSearch(ForceFieldManager& FFM, double MAXDIST,
                                            double LAMBDAMAX, bool *gpu_safestate) {

    //reset safe mode
    _safeMode = false;
    sconvergencecheck = true;
    //calculate first lambda
    double lambda = LAMBDAMAX;
//    std::cout<<lambda<<endl;
    std::cout<<"safe 0"<<endl;
#ifndef CUDAACCL //SERIAL
    cconvergencecheck = new bool[1];
    cconvergencecheck[0] = true;
//#else
//    sconvergencecheck = false;
#endif
//prepare for ping pong optimization
#ifdef CUDAACCL
    nvtxRangePushA("Evcreatesafe");
    CUDAcommon::handleerror(cudaStreamCreate(&s1));
    CUDAcommon::handleerror(cudaStreamCreate(&s2));
    CUDAcommon::handleerror(cudaStreamCreate(&s3));
    CUDAcommon::handleerror(cudaEventCreate(&e1));
    CUDAcommon::handleerror(cudaEventCreate(&e2));
    nvtxRangePop();
    h_stop[0] = false;
    cconvergencecheck = h_stop;
    sp1 = &s1;
    sp2 = &s2;
    ep1 = &e1;
    ep2 = &e2;
    g_s1 = g_stop1;
    g_s2 = g_stop2;
    //prep for safe backtracking.
    if(gpu_params ==NULL){
        double params[5];
        params[0] = BACKTRACKSLOPE;
        params[1] = LAMBDAREDUCE;
        params[2] = LAMBDATOL;
        params[3] = LAMBDAMAX;
        params[4] = MAXDIST;
        CUDAcommon::handleerror(cudaMalloc((void **) &gpu_params, 5 * sizeof(double)));
        CUDAcommon::handleerror(cudaMemcpy(gpu_params, params, 5 * sizeof(double),
                                           cudaMemcpyHostToDevice));
    }
    //initialize lambda search
    CUDAinitializeLambda(*sp1, g_s1, g_s2, gpu_safestate, gpu_state);
    CUDAcommon::handleerror(cudaEventRecord(*ep1, *sp1));
    auto cvars = CUDAcommon::getCUDAvars();
    cvars.streamvec.clear();
    CUDAcommon::cudavars = cvars;
    //check if converged.
    nvtxRangePushA("safeConv");
    CUDAcommon::handleerror(cudaStreamWaitEvent(s3, *ep1, 0));
    CUDAcommon::handleerror(cudaMemcpyAsync(h_stop, g_s2, sizeof(bool), cudaMemcpyDeviceToHost, s3));
    nvtxRangePop();
//        double cudalambda[1];
//    CUDAcommon::handleerror(cudaMemcpy(cudalambda, CUDAcommon::getCUDAvars().gpu_lambda, sizeof(double),
//                                       cudaMemcpyDeviceToHost));
//    std::cout<<lambda<<" "<<cudalambda[0]<<" "<<cconvergencecheck[0]<< sconvergencecheck<<endl;
//    std::cout<<endl;
#endif
    nvtxRangePushA("Energy S 0");
    double currentEnergy = FFM.computeEnergy(coord, force, 0.0);
    nvtxRangePop();
#ifdef CUDAACCL
    //wait for energies to be calculated
    nvtxRangePushA("Safebacksync");
//    std::cout<<"Total energy streams "<<CUDAcommon::getCUDAvars().streamvec.size()<<endl;
    for(auto strm:CUDAcommon::getCUDAvars().streamvec)
    {CUDAcommon::handleerror(cudaStreamSynchronize(*strm));}
    nvtxRangePop();
    setcurrentenergy<<<1,1,0,*sp1>>>(CUDAcommon::getCUDAvars().gpu_energy, g_currentenergy, CUDAcommon::getCUDAvars()
            .gpu_lambda, gpu_initlambdalocal);
    CUDAcommon::handleerror(cudaGetLastError(),"setcurrentenergy", "CGMethod.cu");
    CUDAcommon::handleerror(cudaStreamSynchronize (*sp1));
    h_stop[0] = false;
    cconvergencecheck = h_stop;
#endif
int iter =0;
    //safe backtracking loop
    while(!(cconvergencecheck[0])||!(sconvergencecheck)) {
    iter++;
#ifdef CUDAACCL
        CUDAcommon::handleerror(cudaStreamWaitEvent(*sp2, *ep1, 0));
        //ping pong swap
        sps = sp1;
        sp1 = sp2;
        sp2 = sps;
        eps = ep1;
        ep1 = ep2;
        ep2 = eps;
//        g_bs = g_b1;
//        g_b1 = g_b2;
//        g_b2 = g_bs;
        g_ss = g_s1;
        g_s1 = g_s2;
        g_s2 = g_ss;
        auto cvars = CUDAcommon::getCUDAvars();
        cvars.streamvec.clear();
        CUDAcommon::cudavars = cvars;
#endif

        //new energy when moved by lambda
        std::cout<<"safe z"<<endl;
        nvtxRangePushA("Energy S");
        double energyLambda = FFM.computeEnergy(coord, force, lambda);
        nvtxRangePop();

#ifdef CUDAACCL
        //wait for energies to be calculated
        nvtxRangePushA("Safebacksync");
        for(auto strm:CUDAcommon::getCUDAvars().streamvec)
            CUDAcommon::handleerror(cudaStreamSynchronize(*strm),"Stream Synchronize", "CGMethod.cu");
        nvtxRangePop();

        if(!(cconvergencecheck[0])){
//            CUDAfindLambda(*sp1, stream_bt, g_s1, g_s2, gpu_safestate, gpu_state);
            //check if converged.
            nvtxRangePushA("safeConv");
            if(cconvergencecheck[0]  == false){

                CUDAcommon::handleerror(cudaStreamWaitEvent(s3, *ep1, 0));
                cudaMemcpyAsync(h_stop, g_s2, sizeof(bool), cudaMemcpyDeviceToHost, s3);
            }
            nvtxRangePop();
//            CUDAcommon::handleerror(cudaMemcpy(cudalambda, CUDAcommon::getCUDAvars().gpu_lambda, sizeof(double),
//                                               cudaMemcpyDeviceToHost));
//            std::cout<<lambda<<" "<<cudalambda[0]<<" "<<cconvergencecheck[0]<< sconvergencecheck<<endl;
//            std::cout<<endl;
        }
#else
        if(!(sconvergencecheck)){
            double energyChange = energyLambda - currentEnergy;

            //return if ok
            if(energyChange <= 0.0) sconvergencecheck = true;
            else
                //reduce lambda
                lambda *= LAMBDAREDUCE;

            //just shake if we cant find an energy min,
            //so we dont get stuck
            if(lambda <= 0.0 || lambda <= LAMBDATOL) {
                lambda = MAXDIST / maxF();
                sconvergencecheck = true;
            }
        }
#endif
    }

#ifdef CUDAACCL
    nvtxRangePushA("Evsyncsafe");
    CUDAcommon::handleerror(cudaStreamSynchronize(s1));
    CUDAcommon::handleerror(cudaStreamSynchronize(s2));
    CUDAcommon::handleerror(cudaStreamSynchronize(s3));
    nvtxRangePop();
    nvtxRangePushA("EvDESTROYsafe");
    CUDAcommon::handleerror(cudaStreamDestroy(s1));
    CUDAcommon::handleerror(cudaStreamDestroy(s2));
    CUDAcommon::handleerror(cudaStreamDestroy(s3));
    CUDAcommon::handleerror(cudaEventDestroy(e1));
    CUDAcommon::handleerror(cudaEventDestroy(e2));
    nvtxRangePop();
#else
    delete cconvergencecheck;
#endif
    std::cout<<"lambda determined in "<<iter<<endl;
    if(cconvergencecheck[0]||sconvergencecheck)
        return lambda;;
}
