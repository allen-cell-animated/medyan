
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

#include "CGMethod.h"

#include <algorithm> // max

#include "ForceFieldManager.h"

#include "CGMethodCUDA.h"
#include "MathFunctions.h"
#include "Structure/Bead.h"
#ifdef CUDAACCL
#ifdef __CUDACC__
#define CUDA_HOSTDEV __host__ __device__
#else
#define CUDA_HOSTDEV
#include "nvToolsExt.h"
#endif
#include <cuda.h>
#include <cuda_runtime.h>
#include "CUDAcommon.h"
#endif
#define ARRAY_SIZE 128
//
#include <vector>
#include <cmath>
#include <ctime>
#include "Bead.h"
#include <ctime>
#include <cstdlib>
#include "cross_check.h"
//
long CGMethod::N = 0;
long CGMethod::Ncyl = 0;
#ifdef CUDAACCL
void CGMethod::CUDAresetlambda(cudaStream_t stream) {
    resetlambdaCUDA<<<1,1,0, stream>>>(CUDAcommon::getCUDAvars().gpu_lambda);
            CUDAcommon::handleerror(cudaGetLastError(), "resetlambdaCUDA", "CGMethod.cu");
}
void CGMethod::CUDAinitializeLambda(cudaStream_t stream, bool *check_in, bool *check_out, bool *Polaksafestate, int
                                    *gpu_state){

//    cudaStream_t  s;
//    CUDAcommon::handleerror(cudaStreamCreate(&s));
//    cudaEvent_t  e;
//    CUDAcommon::handleerror(cudaEventCreate(&e));


//    maxFCUDA<<<1,1, 0, s>>>(CUDAcommon::getCUDAvars().gpu_forceAux, gpu_nint, gpu_fmax);
//    cudaStreamSynchronize(s);

//    maxFCUDAred<<<1,3, 3*sizeof(double), s>>>(CUDAcommon::getCUDAvars().gpu_forceAux, gpu_nint, gpu_fmax);
//    cudaStreamSynchronize(s);


//    CUDAcommon::handleerror(cudaDeviceSynchronize());
//    std::cout<<"======"<<endl;
//    CUDAcommon::handleerror(cudaEventRecord(e,s));
//    CUDAcommon::handleerror(cudaGetLastError(), "maxFCUDA", "CGMethod.cu");

////    CUDAcommon::handleerror(cudaStreamWaitEvent(stream,e,0));

////    CUDAcommon::handleerror(cudaEventDestroy(e));

//    CUDAcommon::handleerror(cudaStreamDestroy(s));

//    auto gpu_lambda = CUDAcommon::getCUDAvars().gpu_lambda;
    auto gpu_energy = CUDAcommon::getCUDAvars().gpu_energy;
    initializeLambdaCUDA<<<1,1,0, stream>>>(check_in, check_out, g_currentenergy, gpu_energy, gpu_initlambdalocal, gpu_fmax,
            gpu_params, Polaksafestate, gpu_state);

//    CUDAcommon::handleerror(cudaStreamSynchronize (stream));

    CUDAcommon::handleerror(cudaGetLastError(), "initializeLambdaCUDA", "CGMethod.cu");
}

//void CGMethod::getmaxFCUDA(double *gpu_forceAux, int *gpu_nint, double *gpu_fmax) {
//    maxFCUDA<<<1,1>>>(CUDAcommon::getCUDAvars().gpu_forceAux, gpu_nint, gpu_fmax);
//    CUDAcommon::handleerror(cudaGetLastError(), "getmaxFCUDA", "CGMethod.cu");
//}
void CGMethod::CUDAfindLambda(cudaStream_t  stream1, cudaStream_t stream2, cudaEvent_t  event, bool *checkin, bool
        *checkout, bool *gpu_safestate, int *gpu_state) {
//ToDo remove stream2 from the list of args.
    auto gpu_energy = CUDAcommon::getCUDAvars().gpu_energy;
    auto gpu_lambda = CUDAcommon::getCUDAvars().gpu_lambda;
    findLambdaCUDA << < 1, 1, 0, stream1 >> > (gpu_energy, g_currentenergy, gpu_FDotFA, gpu_fmax, gpu_lambda,
            gpu_params, checkin, checkout, gpu_safestate, gpu_state);
    CUDAcommon::handleerror(cudaEventRecord(event, stream1));
#ifdef ALLSYNC
    cudaDeviceSynchronize();
#endif
/*    findLambdaCUDA2 << < 1, 1, 0, stream1 >> > (gpu_fmax, gpu_lambda, gpu_params,
            checkin, checkout, gpu_safestate,
            gpu_state);
    CUDAcommon::handleerror(cudaGetLastError(), "findLambdaCUDA", "CGMethod.cu")*/;
}
//void CGMethod::CUDAprepforbacktracking(cudaStream_t stream, bool *check_in, bool *check_out){

//    cudaStream_t  s;
//    CUDAcommon::handleerror(cudaStreamCreate(&s));
//    cudaEvent_t  e;
//    CUDAcommon::handleerror(cudaEventCreate(&e));

//    maxFCUDA<<<1,1, 0, s>>>(CUDAcommon::getCUDAvars().gpu_forceAux, gpu_nint, gpu_fmax);
//    CUDAcommon::handleerror(cudaEventRecord(e,s));
//    CUDAcommon::handleerror(cudaGetLastError());


//    CUDAcommon::handleerror(cudaStreamWaitEvent(stream,e,0));


//    CUDAcommon::handleerror(cudaEventDestroy(e));
//    CUDAcommon::handleerror(cudaStreamDestroy(s));

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

    allFADotFCUDA<<<blocksnthreads[0], blocksnthreads[1],0,stream>>>(CUDAcommon::getCUDAvars().gpu_force,
            CUDAcommon::getCUDAvars().gpu_force ,gpu_g, gpu_nint);
    CUDAcommon::handleerror(cudaGetLastError(), "allFADotFCUDA", "CGMethod.cu");
//    addvector<<<1,1,0,stream>>>(gpu_g, gpu_nint, gpu_FDotF);
//    cudaStreamSynchronize(stream);
//    addvectorred<<<1,200,200 * sizeof(double),stream>>>(gpu_g, gpu_nint, gpu_FDotF);
//    double Sum[1];
//        CUDAcommon::handleerror(cudaMemcpy(Sum, gpu_FDotF, sizeof(double), cudaMemcpyDeviceToHost));
    resetdoublevariableCUDA<<<1,1,0,stream>>>(gpu_FDotF);
    addvectorredcgm<<<bntaddvector.at(2),bntaddvector.at(3), bntaddvector.at(3) * sizeof(double),stream>>>(gpu_g,
            gpu_nint, gpu_FDotF);
//    double Sum2[1];
//    CUDAcommon::handleerror(cudaMemcpy(Sum2, gpu_FDotF, sizeof(double), cudaMemcpyDeviceToHost));
//    std::cout<<Sum[0]<<" "<<Sum2[0]<<endl;
//    cudaStreamSynchronize(stream);
    CUDAcommon::handleerror(cudaGetLastError(), "allFADotFCUDA", "CGMethod.cu");

}
void CGMethod::CUDAallFADotFA(cudaStream_t stream){

    allFADotFCUDA<<<blocksnthreads[0], blocksnthreads[1],0,stream>>>(CUDAcommon::getCUDAvars().gpu_forceAux,
            CUDAcommon::getCUDAvars().gpu_forceAux ,gpu_g, gpu_nint);
    CUDAcommon::handleerror(cudaGetLastError(), "allFADotFCUDA", "CGMethod.cu");
//    addvector<<<1,1,0,stream>>>(gpu_g, gpu_nint, gpu_FADotFA);
//    cudaStreamSynchronize(stream);
//    addvectorred<<<1,200,200* sizeof(double),stream>>>(gpu_g, gpu_nint, gpu_FADotFA);
//    cudaStreamSynchronize(stream);
    resetdoublevariableCUDA<<<1,1,0,stream>>>(gpu_FADotFA);
    addvectorredcgm<<<bntaddvector.at(2),bntaddvector.at(3), bntaddvector.at(3) * sizeof(double),stream>>>(gpu_g,
            gpu_nint, gpu_FADotFA);
    CUDAcommon::handleerror(cudaGetLastError(), "allFADotFCUDA", "CGMethod.cu");

}
void CGMethod::CUDAallFADotFAP(cudaStream_t stream){

    allFADotFCUDA<<<blocksnthreads[0], blocksnthreads[1],0,stream>>>(CUDAcommon::getCUDAvars().gpu_forceAux,
            CUDAcommon::getCUDAvars().gpu_forceAuxP ,gpu_g, gpu_nint);
    CUDAcommon::handleerror(cudaGetLastError(), "allFADotFCUDA", "CGMethod.cu");
//    addvector<<<1,1,0,stream>>>(gpu_g, gpu_nint, gpu_FADotFAP);
//    cudaStreamSynchronize(stream);
//    addvectorred<<<1,200,200 * sizeof(double),stream>>>(gpu_g, gpu_nint, gpu_FADotFAP);
    resetdoublevariableCUDA<<<1,1,0,stream>>>(gpu_FADotFAP);
    addvectorredcgm<<<bntaddvector.at(2),bntaddvector.at(3), bntaddvector.at(3) * sizeof(double),stream>>>(gpu_g,
            gpu_nint, gpu_FADotFAP);
//    cudaStreamSynchronize(stream);
    CUDAcommon::handleerror(cudaGetLastError(), "allFADotFCUDA", "CGMethod.cu");

}
void CGMethod::CUDAallFDotFA(cudaStream_t stream){

    allFADotFCUDA<<<blocksnthreads[0], blocksnthreads[1],0,stream>>>(CUDAcommon::getCUDAvars().gpu_force,
            CUDAcommon::getCUDAvars().gpu_forceAux ,gpu_g, gpu_nint);
    CUDAcommon::handleerror(cudaGetLastError(), "allFADotFCUDA", "CGMethod.cu");
//    addvector<<<1,1,0,stream>>>(gpu_g, gpu_nint, gpu_FDotFA);
//    cudaStreamSynchronize(stream);
//    addvectorred<<<1,200,200* sizeof(double),stream>>>(gpu_g, gpu_nint, gpu_FDotFA);
//    cudaStreamSynchronize(stream);
    resetdoublevariableCUDA<<<1,1,0,stream>>>(gpu_FDotFA);
    addvectorredcgm<<<bntaddvector.at(2),bntaddvector.at(3), bntaddvector.at(3) * sizeof(double),stream>>>(gpu_g,
            gpu_nint, gpu_FDotFA);
    CUDAcommon::handleerror(cudaGetLastError(), "allFADotFCUDA", "CGMethod.cu");

}

void CGMethod::CUDAshiftGradient(cudaStream_t stream, bool *Mcheckin) {
    shiftGradientCUDA<<<blocksnthreads[0], blocksnthreads[1],0, stream>>>(CUDAcommon::getCUDAvars().gpu_force,
            CUDAcommon::getCUDAvars().gpu_forceAux, gpu_nint, gpu_FADotFA, gpu_FADotFAP, gpu_FDotF, Mcheckin);
}

void CGMethod::CUDAshiftGradientifSafe(cudaStream_t stream, bool *Mcheckin, bool *Scheckin){
    shiftGradientCUDAifsafe<<<blocksnthreads[0], blocksnthreads[1],0, stream>>>(CUDAcommon::getCUDAvars().gpu_force, CUDAcommon::getCUDAvars().gpu_forceAux, gpu_nint,
                            Mcheckin, Scheckin);
    CUDAcommon::handleerror(cudaGetLastError(),"CUDAshiftGradientifSafe", "CGMethod.cu");
}

//void CGMethod::CUDAgetPolakvars(bool calc_safestate,cudaStream_t streamcalc, double* gpu_GRADTOL, bool *gminstatein,
//                                    bool *gminstateout, bool *gsafestateout, volatile bool *cminstate){
////    state[0] = false;
////    state[1] = false;
//    if(cminstate[0] == true) {

////        maxFCUDA << < 1, 1, 0, streamcalc >> > (CUDAcommon::getCUDAvars().gpu_forceAux, gpu_nint, gpu_fmax);
//        maxFCUDAred<<<1,3, 3*sizeof(double), streamcalc>>>(CUDAcommon::getCUDAvars().gpu_forceAux, gpu_nint, gpu_fmax);
////        CUDAcommon::handleerror(cudaDeviceSynchronize());
////        std::cout<<"======"<<endl;
//        CUDAcommon::handleerror(cudaGetLastError(), "maxFCUDA", "CGMethod.cu");

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

//        maxFCUDA << < 1, 1, 0, streamcalc >> > (CUDAcommon::getCUDAvars().gpu_forceAux, gpu_nint, gpu_fmax);
//        maxFCUDAred<<<1,3, 3*sizeof(double), streamcalc>>>(CUDAcommon::getCUDAvars().gpu_forceAux, gpu_nint, gpu_fmax);
//        cudaStreamSynchronize(streamcalc);

        //@{ V2
//        allFADotFCUDA<<<blocksnthreads[0], blocksnthreads[1],0,streamcalc>>>(CUDAcommon::getCUDAvars().gpu_forceAux,
//                CUDAcommon::getCUDAvars().gpu_forceAux ,gpu_maxF, gpu_nint);
//        CUDAcommon::handleerror(cudaGetLastError(), "allFADotFACUDA", "CGMethod.cu");
//        maxFCUDAredv2<<<1,512,512*sizeof(double), streamcalc>>>(gpu_maxF, gpu_nint,
//                gpu_fmax);
        //@}
        //Test
//        CUDAcommon::handleerror(cudaDeviceSynchronize());
//        double maxFv1[1];
//        cudaMemcpy(maxFv1, gpu_fmax,  sizeof(double), cudaMemcpyDeviceToHost);
//        std::cout<<"v1 maxF "<<maxFv1[0]<<endl;
//        double *gpu_fmax2;
//        CUDAcommon::handleerror(cudaMalloc((void **)&gpu_fmax2, sizeof(double)));
/*#ifdef CUDATIMETRACK
        cudaStream_t  streamcalc2;
        cudaStreamCreate(&streamcalc2);
        streamcalc = streamcalc2;
        chrono::high_resolution_clock::time_point tbegin, tend;
        tbegin = chrono::high_resolution_clock::now();
#endif*/
        //@{ V3
        //TODO combine with FADOTFA calculation by making it write it to gpu_maxF before
        // adding.
        allFADotFCUDA<<<blocksnthreads[0], blocksnthreads[1],0,streamcalc>>>(CUDAcommon::getCUDAvars().gpu_forceAux,
                CUDAcommon::getCUDAvars().gpu_forceAux ,gpu_maxF, gpu_nint);
        resetdoublevariableCUDA<<<1,1,0,streamcalc>>>(gpu_fmax);
        resetintvariableCUDA<<<1,1,0,streamcalc>>>(gpu_mutexlock);
        maxFCUDAredv3<<<bntaddvector.at(2),bntaddvector.at(3), bntaddvector.at(3) *
                sizeof(double),streamcalc>>>(gpu_maxF, gpu_nint, gpu_fmax, gpu_mutexlock);
        CUDAcommon::handleerror(cudaGetLastError(), "maxFCUDA", "CGMethod.cu");
        getminimizestateCUDA << < 1, 1, 0, streamcalc >> > (gpu_fmax, gpu_GRADTOL, gminstatein, gminstateout);
        CUDAcommon::handleerror(cudaGetLastError(), "getminimizestateCUDA", "CGMethod.cu");
        //@}
/*#ifdef CUDATIMETRACK
        cudaStreamSynchronize(streamcalc);
        tend = chrono::high_resolution_clock::now();
        chrono::duration<double> elapsed_runs1(tend - tbegin);
        std::cout<<"CUDA maxF "<<elapsed_runs1.count()<<endl;
#endif
#ifdef CUDATIMETRACK
        tbegin = chrono::high_resolution_clock::now();
#endif
        maxF();
        auto x = maxF()>maxF();
#ifdef CUDATIMETRACK
        cudaStreamSynchronize(streamcalc);
        tend = chrono::high_resolution_clock::now();
        chrono::duration<double> elapsed_runs2(tend - tbegin);
        std::cout<<"SERL maxF "<<elapsed_runs2.count()<<endl;
#endif*/


//        CUDAcommon::handleerror(cudaDeviceSynchronize());
//        cudaMemcpy(maxFv1, gpu_fmax2,  sizeof(double), cudaMemcpyDeviceToHost);
//        std::cout<<"v2 maxF "<<maxFv1[0]<<endl;
//        cudaFree(gpu_fmax2);
        //Test ends

//                cout<<"MaxF algorithm is not accurate. Redo algorithm. Exiting"<<endl;
//                exit(EXIT_FAILURE);

//        cudaStreamSynchronize(streamcalc);
//        CUDAcommon::handleerror(cudaDeviceSynchronize());
//        std::cout<<"======"<<endl;
/*        CUDAcommon::handleerror(cudaGetLastError(), "maxFCUDA", "CGMethod.cu");
        getminimizestateCUDA << < 1, 1, 0, streamcalc >> > (gpu_fmax, gpu_GRADTOL, gminstatein, gminstateout);
        CUDAcommon::handleerror(cudaGetLastError(), "getminimizestateCUDA", "CGMethod.cu");*/
//        CUDAcommon::handleerror(cudaStreamSynchronize(streamcalc));
    }
//    if(calc_safestate){
//        CUDAallFDotFA(streamcalc);
//        getsafestateCUDA<<<1,1,0,streamcalc>>>(gpu_FDotFA, gpu_FDotF, gpu_FADotFA, gsafestateout);
//        CUDAcommon::handleerror(cudaGetLastError(), "getsafestateCUDA", "CGMethod.cu");
//    }
    CUDAcommon::handleerror(cudaGetLastError(),"CUDAgetPolakvars", "CGMethod.cu");
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

    moveBeadsCUDA<<<blocksnthreads[0], blocksnthreads[1],0, stream>>>(gpu_coord, gpu_force, gpu_lambda, gpu_nint,
            gpu_checkin);

    CUDAcommon::handleerror(cudaGetLastError(),"moveBeadsCUDA", "CGMethod.cu");
}

void CGMethod::CUDAinitializePolak(cudaStream_t stream, bool *minstatein, bool *minstateout, bool *safestatein, bool
        *safestateout){
    CUDAallFDotFA(stream);
    initializePolak<<<1,1,0,stream>>>(minstatein, minstateout, safestatein, safestateout);
    CUDAcommon::handleerror(cudaGetLastError(),"CUDAinitializePolak", "CGPolakRibiereMethod.cu");
}

//double CGMethod::gpuFDotF(double *f1,double *f2){
//
//    allFADotFCUDA<<<blocksnthreads[0], blocksnthreads[1]>>>(f1, f2 ,gpu_g, gpu_nint);
//    CUDAcommon::handleerror(cudaGetLastError(),"allFADotFCUDA", "CGMethod.cu");
////    addvector<<<1,1>>>(gpu_g, gpu_nint, gSum);
//    addvectorred<<<1,200,200* sizeof(double)>>>(gpu_g, gpu_nint, gSum);
//    CUDAcommon::handleerror(cudaGetLastError(),"allFADotFCUDA", "CGMethod.cu");
//
////    CUDAcommon::handleerror( cudaPeekAtLastError() );
////    CUDAcommon::handleerror(cudaDeviceSynchronize());
//
//    double g[1];
//    CUDAcommon::handleerror(cudaMemcpy(g, gSum, sizeof(double),
//                                       cudaMemcpyDeviceToHost));
//
//
////    double g[N/3];
////    CUDAcommon::handleerror(cudaMemcpy(g, gpu_g, N/3 * sizeof(double),
////                                       cudaMemcpyDeviceToHost));
////    CUDAcommon::handleerror(cudaFree(gpu_g));
////    double sum=0.0;
////    for(auto i=0;i<N/3;i++)
////        sum+=g[i];
//    return g[0];
//}
#endif
double CGMethod::allFDotF()
{
    return mathfunc::dot(Bead::getDbData().forces, Bead::getDbData().forces);
}

double CGMethod::allFADotFA()
{
//#ifdef CUDAACCL

//    auto g_cuda = gpuFDotF(CUDAcommon::getCUDAvars().gpu_forceAux,CUDAcommon::getCUDAvars().gpu_forceAux);

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
    return mathfunc::dot(Bead::getDbData().forcesAux, Bead::getDbData().forcesAux);
}

double CGMethod::allFADotFAP()
{
    return mathfunc::dot(Bead::getDbData().forcesAux, Bead::getDbData().forcesAuxP);
}

double CGMethod::allFDotFA()
{
//#ifdef CUDAACCL
//    auto g_cuda = gpuFDotF(CUDAcommon::getCUDAvars().gpu_force,CUDAcommon::getCUDAvars().gpu_forceAux);
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
    return mathfunc::dot(Bead::getDbData().forces, Bead::getDbData().forcesAux);
}

double CGMethod::maxF() {
    double mag2Max = 0.0;
    for(auto x : Bead::getDbData().forcesAux) {
        mag2Max = std::max(mag2Max, mathfunc::magnitude2(x));
    }

    return std::sqrt(mag2Max);
}

Bead* CGMethod::maxBead() {

    double maxF2 = 0.0;
    double currentF2;
    long index = 0;
#ifdef SERIAL
    const std::size_t numBeads = Bead::numBeads();
    for (size_t i = 0; i < numBeads; ++i) {
        currentF2 = mathfunc::magnitude2(Bead::getDbData().forcesAux[i]);
        if(currentF2 > maxF2) {
            index = i;
            maxF2 = currentF2;
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
            index = (i - i%3)/3;
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

//    std::cout<<"3N "<<N<<endl;
    const std::size_t num = Bead::getDbData().coords.size_raw();
    for(size_t i = 0; i < num; ++i)
        Bead::getDbData().coords.value[i] += d * Bead::getDbData().forces.value[i];
}

void CGMethod::shiftGradient(double d)
{
    const std::size_t num = Bead::getDbData().coords.size_raw();
    for (size_t i = 0; i < num; i ++)
        Bead::getDbData().forces.value[i]
            = Bead::getDbData().forcesAux.value[i] + d * Bead::getDbData().forces.value[i];
}

void CGMethod::printForces()
{
    cout << "Print Forces" << endl;
    for(auto b: Bead::getBeads()) {

        for (int i = 0; i<3; i++)
            cout << b->coordinate()[i] << "  "<<
                 b->force()[i] <<"  "<<b->forceAux()[i]<<endl;
    }
    cout << "End of Print Forces" << endl;
}

void CGMethod::startMinimization() {
    //COPY BEAD DATA
#ifdef CUDATIMETRACK
    chrono::high_resolution_clock::time_point tbegin, tend;
    tbegin = chrono::high_resolution_clock::now();
#endif
    Ncyl = Cylinder::getCylinders().size();

#ifdef CUDATIMETRACK
    tend= chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed_runst(tend - tbegin);
    CUDAcommon::cudatime.Tstartmin = elapsed_runst.count();
    std::cout<<"Start conv to vec time taken (s) "<<elapsed_runst.count()<<endl;
#endif
#ifdef CUDAACCL
#ifdef CUDATIMETRACK
    tbegin = chrono::high_resolution_clock::now();
#endif
    //Start stream
    if(stream_startmin == NULL || !(CUDAcommon::getCUDAvars().conservestreams))
        CUDAcommon::handleerror(cudaStreamCreate(&stream_startmin));
    int nDevices;
//    cudaDeviceProp prop;
    cudaGetDeviceCount(&nDevices);
    if(nDevices>1){
        cout<<"Code not configured for multiple devices. Exiting..."<<endl;
        exit(EXIT_FAILURE);
    }

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
    cylinder* gpu_cylindervec;
    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_cylindervec, Ncyl*sizeof(cylinder)));

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
    CUDAcommon::handleerror(cudaMalloc(&g_stop1, sizeof(bool)));
    CUDAcommon::handleerror(cudaMalloc(&g_stop2, sizeof(bool)));
    CUDAcommon::handleerror(cudaHostAlloc(&h_stop, sizeof(bool), cudaHostAllocDefault));

    //@
    //Store the pointers so they can be tracked while calculating energies.
    CUDAcommon::cudavars.backtrackbools.clear();
    CUDAcommon::cudavars.backtrackbools.push_back(g_stop1);
    CUDAcommon::cudavars.backtrackbools.push_back(g_stop2);

//    CUDAcommon::handleerror(cudaHostAlloc((void**)&convergencecheck, 3 * sizeof(bool), cudaHostAllocDefault));
//    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_convergencecheck, 3 * sizeof(bool)));

//    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_lambda, sizeof(double))); REPEAT.
    CUDAcommon::handleerror(cudaMemcpyAsync(gpu_coord, coord, N*sizeof(double),
                                        cudaMemcpyHostToDevice, stream_startmin));
    CUDAcommon::handleerror(cudaMemcpyAsync(gpu_force, f, N*sizeof(double),
                                        cudaMemcpyHostToDevice, stream_startmin));
    CUDAcommon::handleerror(cudaMemcpyAsync(gpu_forceAux, f, N*sizeof(double),
                                        cudaMemcpyHostToDevice, stream_startmin));
    CUDAcommon::handleerror(cudaMemcpyAsync(gpu_forceAuxP, f, N*sizeof(double),
                                        cudaMemcpyHostToDevice, stream_startmin));
    bool dummy[1];dummy[0] = true;
    CUDAcommon::handleerror(cudaMemcpyAsync(gpu_btstate, dummy, sizeof(bool),
                                        cudaMemcpyHostToDevice, stream_startmin));

    CUDAcommon::handleerror(cudaMemcpyAsync(gpu_cylindervec, cylindervec, Ncyl*sizeof
                                                    (cylinder),
                                            cudaMemcpyHostToDevice, stream_startmin));
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
    cvars.gpu_cylindervec = gpu_cylindervec;
    CUDAcommon::cudavars=cvars;
//SET CERTAIN GPU PARAMETERS SET FOR EASY ACCESS DURING MINIMIZATION._
//    int THREADSPERBLOCK;
//    cudaDeviceProp prop;
//    cudaGetDeviceProperties(&prop, 0);
//    THREADSPERBLOCK = prop.maxThreadsPerBlock;
    //@{ Reduction Add variables
    bntaddvector.clear();
    bntaddvector = getaddred2bnt(N/3);
    int M = bntaddvector.at(0);
    vector<double> zerovec(M);
    fill(zerovec.begin(),zerovec.begin()+M,0.0);
    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_g, M * sizeof(double)));
    CUDAcommon::handleerror(cudaMemcpyAsync(gpu_g, zerovec.data(),
                           M * sizeof(double), cudaMemcpyHostToDevice, stream_startmin));
    /*CUDAcommon::handleerror(cudaMemsetAsync(gpu_g, 0, M * sizeof(double), stream_startmin));*/
    //MaxF
    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_maxF, M * sizeof(double)));
    CUDAcommon::handleerror(cudaMemcpyAsync(gpu_maxF, zerovec.data(),
                                            M * sizeof(double), cudaMemcpyHostToDevice, stream_startmin));
    /*CUDAcommon::handleerror(cudaMemsetAsync(gpu_maxF, 0, M * sizeof(double), stream_startmin));*/
    int THREADSPERBLOCK = bntaddvector.at(1);
    //@}

    int nint[1]; nint[0]=CGMethod::N/3;
    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_nint, sizeof(int)));
    CUDAcommon::handleerror(cudaMemcpyAsync(gpu_nint, nint, sizeof(int),
                                        cudaMemcpyHostToDevice, stream_startmin));
    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_state, sizeof(int)));
    blocksnthreads.push_back(CGMethod::N/(3*THREADSPERBLOCK) + 1);
    if(blocksnthreads[0]==1) blocksnthreads.push_back(CGMethod::N/3);
    else blocksnthreads.push_back(THREADSPERBLOCK);
    auto maxthreads = 8 * THREADSPERBLOCK;

    //@{maxFredv3
    int state[1];state[0] = 0;
    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_mutexlock, sizeof(int)));
    CUDAcommon::handleerror(cudaMemcpyAsync(gpu_mutexlock, state, sizeof(int),
                                        cudaMemcpyHostToDevice, stream_startmin));
    //Synchronize
    CUDAcommon::handleerror(cudaStreamSynchronize(stream_startmin),"CGMethod.cu",
                            "startMinimization");

#ifdef CUDATIMETRACK
    tend= chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed_run(tend - tbegin);
    CUDAcommon::cudatime.Tstartmin = elapsed_run.count();
    std::cout<<"start min time taken (s) "<<elapsed_run.count()<<endl;
#endif

#ifdef CUDATIMETRACK
    CUDAcommon::cudatime.Tlambdap.clear();
    CUDAcommon::cudatime.Tlambdapcount.clear();
    CUDAcommon::cudatime.Tlambdap.push_back(0);
    CUDAcommon::cudatime.Tlambdap.push_back(0);
    CUDAcommon::cudatime.Tlambdap.push_back(0);
    CUDAcommon::cudatime.Tlambdapcount.push_back(0);
    CUDAcommon::cudatime.Tlambdapcount.push_back(0);
    CUDAcommon::cudatime.Tlambdapcount.push_back(0);
    //
    CUDAcommon::serltime.Tlambdap.clear();
    CUDAcommon::serltime.Tlambdapcount.clear();
    CUDAcommon::serltime.Tlambdap.push_back(0);
    CUDAcommon::serltime.Tlambdap.push_back(0);
    CUDAcommon::serltime.Tlambdap.push_back(0);
    CUDAcommon::serltime.Tlambdapcount.push_back(0);
    CUDAcommon::serltime.Tlambdapcount.push_back(0);
    CUDAcommon::serltime.Tlambdapcount.push_back(0);
#endif
    //@}
    //addvectorred2@{

//    int blocks, threads;
//    if(M > THREADSPERBLOCK){
//        if(M > maxthreads) {
//            blocks = 8;
//            threads = THREADSPERBLOCK;
//        }
//        else if(M > THREADSPERBLOCK){
//            blocks = M /(4 * THREADSPERBLOCK) +1;
//            threads = THREADSPERBLOCK;
//        }
//    }
//    else
//    { blocks = 1; threads = M/4;}
//    std::cout<<blocks<<" "<<threads<<" "<<M<<" "<<N/3<<" "<<maxthreads<<" "<<THREADSPERBLOCK<<endl;
//    bntaddvector.clear();
//    bntaddvector.push_back(blocks);
//    bntaddvector.push_back(threads);
//    CUDAcommon::handleerror(cudaMalloc((void **) &gSum, sizeof(double)));
//    CUDAcommon::handleerror(cudaMalloc((void **) &gSum2, sizeof(double)));
    //@}
//    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_g, N/3 * sizeof(double)));

    //Memory alloted
    //@{
//    size_t allocmem = 0;
//    allocmem += (4*N + 9 + M)*sizeof(double) + 6 * sizeof(bool) + 6 * sizeof(int) + 200 * sizeof(char);
//    auto c = CUDAcommon::getCUDAvars();
//    c.memincuda += allocmem;
//    CUDAcommon::cudavars = c;
//    std::cout<<"Total allocated memory KB"<<c.memincuda/1024<<endl;
//    std::cout<<"Memory allocated "<< allocmem/1024<<"Memory freed 0"<<endl;
    //@}


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
#ifdef CUDATIMETRACK
    chrono::high_resolution_clock::time_point tbegin, tend;
    tbegin = chrono::high_resolution_clock::now();
#endif
#ifdef CUDAACCL

    CUDAcommon::handleerror(cudaMemcpy(coord, CUDAcommon::getCUDAvars().gpu_coord, N *
                            sizeof(double), cudaMemcpyDeviceToHost));
    CUDAcommon::handleerror(cudaMemcpy(force, CUDAcommon::getCUDAvars().gpu_force, N *
                            sizeof(double), cudaMemcpyDeviceToHost));
//    CUDAcommon::handleerror(cudaMemcpy(forceAux, CUDAcommon::getCUDAvars().gpu_forceAux, N *
//                            sizeof(double), cudaMemcpyDeviceToHost));

    #endif

//    deallocate();
#ifdef CUDAACCL
    bool deletecndn = true;
#ifdef CUDAACCL_NLS
    deletecndn = false;
#endif
    if(deletecndn) {
        CUDAcommon::handleerror(cudaFree(CUDAcommon::getCUDAvars().gpu_coord));
        CUDAcommon::handleerror(cudaFree(CUDAcommon::getCUDAvars().gpu_cylindervec));
    }
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
//    CUDAcommon::handleerror(cudaFree(gSum));
//    CUDAcommon::handleerror(cudaFree(gSum2));
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
    CUDAcommon::handleerror(cudaFree(gpu_mutexlock));
    blocksnthreads.clear();
    if(!(CUDAcommon::getCUDAvars().conservestreams))
        CUDAcommon::handleerror(cudaStreamDestroy(stream_startmin));

    //TODO cross check later
//    CUDAcommon::handleerror(cudaFree(CUDAcommon::getCUDAvars().motorparams));

//    CUDAcommon::getCUDAvars().gpu_coord = NULL;
//    CUDAcommon::getCUDAvars().gpu_force = NULL;
//    CUDAcommon::getCUDAvars().gpu_forceAux = NULL;
//    CUDAcommon::getCUDAvars().gpu_lambda = NULL;

    //Memory alloted
    //@{
//    size_t allocmem = 0;
//    allocmem += (4*N + 9 +  bntaddvector.at(0))*sizeof(double) + 6 * sizeof(bool) + 6 * sizeof(int) + 200 * sizeof(char);
//    auto c = CUDAcommon::getCUDAvars();
//    c.memincuda -= allocmem;
//    CUDAcommon::cudavars = c;
//    std::cout<<"Total allocated memory "<<c.memincuda/1024<<endl;
//    std::cout<<"Memory allocated 0 . Memory freed "<<allocmem/1024<<endl;
    //@}

//    size_t free, total;
//    CUDAcommon::handleerror(cudaMemGetInfo(&free, &total));
//    fprintf(stdout,"\t### After Min Available VRAM : %g Mo/ %g Mo(total)\n\n",
//            free/1e6, total/1e6);
//
//    cudaFree(0);
//
//    CUDAcommon::handleerror(cudaMemGetInfo(&free, &total));
//    fprintf(stdout,"\t### Available VRAM : %g Mo/ %g Mo(total)\n\n",
//            free/1e6, total/1e6);
#endif
#ifdef CUDATIMETRACK
    tend= chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed_run(tend - tbegin);
    CUDAcommon::cudatime.Tstartmin = elapsed_run.count();
    std::cout<<"end min time taken (s) "<<elapsed_run.count()<<endl;
#endif
}

#ifdef CUDAACCL
double CGMethod::backtrackingLineSearchCUDA(ForceFieldManager& FFM, double MAXDIST,
                                        double LAMBDAMAX, bool *gpu_safestate) {
#ifdef CUDATIMETRACK
    chrono::high_resolution_clock::time_point tbegin, tend;
    CUDAcommon::cudatime.Tlambdapcount.at(0)++;
    tbegin = chrono::high_resolution_clock::now();
#endif
    //@{ Lambda phase 1
    double lambda;
    h_stop[0] = false;
    if(s1 == NULL || !(CUDAcommon::getCUDAvars().conservestreams))
        CUDAcommon::handleerror(cudaStreamCreate(&s1));
    if(s2 == NULL || !(CUDAcommon::getCUDAvars().conservestreams))
        CUDAcommon::handleerror(cudaStreamCreate(&s2));
    if(s3 == NULL || !(CUDAcommon::getCUDAvars().conservestreams))
        CUDAcommon::handleerror(cudaStreamCreate(&s3));
    if(e1 == NULL || !(CUDAcommon::getCUDAvars().conservestreams))
        CUDAcommon::handleerror(cudaEventCreate(&e1));
    if(e2 == NULL || !(CUDAcommon::getCUDAvars().conservestreams))
        CUDAcommon::handleerror(cudaEventCreate(&e2));
    sp1 = &s1;
    sp2 = &s2;
    ep1 = &e1;
    ep2 = &e2;
    g_s1 = g_stop1;
    g_s2 = g_stop2;
    //prep for backtracking.
    if(gpu_params == NULL){
        //TODO move gpu_params copy permanently out of the function.
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
    CUDAresetlambda(*sp1);//set lambda to zero.
    if(e == NULL || !(CUDAcommon::getCUDAvars().conservestreams))  {
        CUDAcommon::handleerror(cudaEventCreate(&e));
    }

    CUDAcommon::handleerror(cudaEventRecord(e, *sp1));
    auto cvars = CUDAcommon::getCUDAvars();
    cvars.streamvec.clear();
    CUDAcommon::cudavars = cvars;
    //initialize lambda search
    CUDAinitializeLambda(*sp1, g_s1, g_s2, gpu_safestate, gpu_state);
    //@} Lambda phase 1
#ifdef CUDATIMETRACK
    tend= chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed_run(tend - tbegin);
    CUDAcommon::cudatime.Tlambdap.at(0) += elapsed_run.count();
#endif
    //Calculate current energy.
    double currentEnergy = FFM.computeEnergy(coord, force, 0.0);
    //wait for energies to be calculated
    for(auto strm:CUDAcommon::getCUDAvars().streamvec)
        CUDAcommon::handleerror(cudaStreamSynchronize(*strm),"backConvSync","CGMethod.cu");
    if(stream_bt == NULL || !(CUDAcommon::getCUDAvars().conservestreams))
        CUDAcommon::handleerror(cudaStreamCreate(&stream_bt),"find lambda", "CGMethod.cu");

#ifdef DETAILEDOUTPUT_ENERGY
//    CUDAcommon::handleerror(cudaDeviceSynchronize());
    double cuda_energy[1];
    CUDAcommon::handleerror(cudaMemcpy(cuda_energy, CUDAcommon::cudavars.gpu_energy,  sizeof(double),
                                       cudaMemcpyDeviceToHost));
    std::cout<<"Total Energy cE pN.nm CUDA "<<cuda_energy[0]<<" SERL "<<currentEnergy<<endl;
    std::cout<<endl;
#endif

#ifdef CUDATIMETRACK
    tbegin = chrono::high_resolution_clock::now();
#endif
    //@{ Lambda phase 1b
    cudaStreamSynchronize(*sp1);
    setcurrentenergy<<<1,1,0,*sp1>>>(CUDAcommon::getCUDAvars().gpu_energy, g_currentenergy, CUDAcommon::getCUDAvars()
            .gpu_lambda, gpu_initlambdalocal);
    CUDAcommon::handleerror(cudaGetLastError(),"setcurrentenergy", "CGMethod.cu");
    cudaStreamSynchronize(*sp1);

    //check if converged.
    //TODO commented coz this line is not needed
//    CUDAcommon::handleerror(cudaStreamWaitEvent(s3, *ep1, 0));
//    CUDAcommon::handleerror(cudaEventRecord(*CUDAcommon::getCUDAvars().event, *sp1));
    CUDAcommon::handleerror(cudaMemcpyAsync(h_stop, g_s2, sizeof(bool), cudaMemcpyDeviceToHost, s3));
//    CUDAcommon::handleerror(cudaStreamSynchronize (*sp1)); CHECK IF NEEDED
    cconvergencecheck = h_stop;
    int iter = 0;
    //@} Lambda phase 1b
#ifdef CUDATIMETRACK
    tend= chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed_run1b(tend - tbegin);
    CUDAcommon::cudatime.Tlambdap.at(0) += elapsed_run1b.count();
#endif

    while(!(cconvergencecheck[0])) {
#ifdef CUDATIMETRACK
        CUDAcommon::cudatime.Tlambdapcount.at(1)++;
        tbegin = chrono::high_resolution_clock::now();
#endif
        //@{ Lambda phase 2
        iter++;
        CUDAcommon::handleerror(cudaStreamWaitEvent(*sp2, *ep1, 0));
        CUDAcommon::handleerror(cudaStreamSynchronize(*sp2));
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
        //@} Lambda phase 2
#ifdef CUDATIMETRACK
        tend= chrono::high_resolution_clock::now();
        chrono::duration<double> elapsed_run2(tend - tbegin);
        CUDAcommon::cudatime.Tlambdap.at(1) += elapsed_run2.count();
#endif

#ifdef SERIAL_CUDACROSSCHECK
        double cuda_lambda[1];
        CUDAcommon::handleerror(cudaDeviceSynchronize(),"CGPolakRibiereMethod.cu","CGPolakRibiereMethod.cu");
        CUDAcommon::handleerror(cudaMemcpy(cuda_lambda, CUDAcommon::cudavars.gpu_lambda,  sizeof(double),
                                           cudaMemcpyDeviceToHost));
        lambda = cuda_lambda[0];
#endif

        //TODO let each forcefield calculate energy IFF conv state = false. That will help
        // them avoid unnecessary iterations.
        //let each forcefield also add energies to two different energy variables.
        double energyLambda = FFM.computeEnergy(coord, force, lambda);

        //wait for energies to be calculated
         for(auto strm:CUDAcommon::getCUDAvars().streamvec) {
            CUDAcommon::handleerror(cudaStreamSynchronize(*strm), "backConvsync", "CGMethod.cu");
        }
#ifdef SERIAL_CUDACROSSCHECK
        for(auto strm:CUDAcommon::getCUDAvars().streamvec) {
            CUDAcommon::handleerror(cudaStreamSynchronize(*strm), "backConvsync", "CGMethod.cu");
        }
        CUDAcommon::handleerror(cudaDeviceSynchronize());
        double cuda_energy[1];
        CUDAcommon::handleerror(cudaMemcpy(cuda_energy, CUDAcommon::cudavars.gpu_energy,  sizeof(double),
                                           cudaMemcpyDeviceToHost));
        std::cout<<"Total Energy EL pN.nm CUDA "<<cuda_energy[0]<<" SERL "
                ""<<energyLambda<<endl;
        std::cout<<endl;
#endif
#ifdef CUDATIMETRACK
        tbegin = chrono::high_resolution_clock::now();
#endif
        //@{ Lambda phase 2
        if(!(cconvergencecheck[0])){
            CUDAcommon::handleerror(cudaStreamSynchronize(stream_bt));
            CUDAfindLambda(*sp1, stream_bt, *ep1, g_s1, g_s2, gpu_safestate, gpu_state);
            CUDAcommon::handleerror(cudaStreamSynchronize(*sp1));
            CUDAcommon::handleerror(cudaStreamSynchronize(stream_bt));
            if(cconvergencecheck[0]  == false){
                CUDAcommon::handleerror(cudaStreamWaitEvent(s3, *ep1, 0));
                CUDAcommon::handleerror(cudaMemcpyAsync(h_stop, g_s2, sizeof(bool), cudaMemcpyDeviceToHost, s3));
            }
        }
        //@Lambda phase 2
#ifdef CUDATIMETRACK
        tend= chrono::high_resolution_clock::now();
        chrono::duration<double> elapsed_run2b(tend - tbegin);
        CUDAcommon::cudatime.Tlambdap.at(1) += elapsed_run2b.count();
#endif
    }
    if(!(CUDAcommon::getCUDAvars().conservestreams))
        CUDAcommon::handleerror(cudaFree(gpu_params), "CudaFree", "CGMethod.cu");
#ifdef CUDATIMETRACK
    CUDAcommon::cudatime.Tlambdapcount.at(2)++;
    tbegin = chrono::high_resolution_clock::now();
#endif
    //@{ Lambda phase 3
    //commented on 18 Sep 2018.
//    correctlambdaCUDA<<<1,1,0, stream_bt>>>(CUDAcommon::getCUDAvars().gpu_lambda, gpu_state, gpu_params);

/*    correctlambdaCUDA<<<1,1,0, *sp1>>>(CUDAcommon::getCUDAvars().gpu_lambda, gpu_state,
            gpu_params);*/

    CUDAcommon::handleerror(cudaStreamSynchronize(stream_bt));
    CUDAcommon::handleerror(cudaStreamSynchronize(s1));
    CUDAcommon::handleerror(cudaStreamSynchronize(s2));
    CUDAcommon::handleerror(cudaStreamSynchronize(s3));
    //@} Lambda phase 3
#ifdef CUDATIMETRACK
    tend= chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed_run3(tend - tbegin);
    CUDAcommon::cudatime.Tlambdap.at(2) += elapsed_run3.count();
#endif
    if(!(CUDAcommon::getCUDAvars().conservestreams))  {
        CUDAcommon::handleerror(cudaStreamDestroy(s1));
        CUDAcommon::handleerror(cudaStreamDestroy(s2));
        CUDAcommon::handleerror(cudaStreamDestroy(s3));
        CUDAcommon::handleerror(cudaStreamDestroy(stream_bt));
        CUDAcommon::handleerror(cudaEventDestroy(e1));
        CUDAcommon::handleerror(cudaEventDestroy(e2));
    }
    std::cout<<"CUDA lambda determined in "<<iter<< " iterations "<<endl;

    if(cconvergencecheck[0]||sconvergencecheck)
        return lambda;

}
#endif // CUDAACCL

double CGMethod::backtrackingLineSearch(ForceFieldManager& FFM, double MAXDIST,
                                        double LAMBDAMAX, bool *gpu_safestate) {

    //@{ Lambda phase 1
    double lambda;
    sconvergencecheck = true;
#ifdef SERIAL //SERIAL
    sconvergencecheck = false;
    cconvergencecheck = new bool[1];
    cconvergencecheck[0] = true;
#endif
#ifdef SERIAL
    double f = maxF();
    //return zero if no forces
    if(f == 0.0){
        lambda = 0.0;
#ifdef DETAILEDOUTPUT_LAMBDA
        std::cout<<"initial_lambda_serial "<<lambda<<endl;
#endif
        sconvergencecheck = true;}
    //calculate first lambda
    lambda = min(LAMBDAMAX, MAXDIST / f);

    //@} Lambda phase 1
#ifdef DETAILEDOUTPUT_LAMBDA
    std::cout<<"SL lambdamax "<<LAMBDAMAX<<" serial_lambda "<<lambda<<" fmax "<<f<<" state "<<sconvergencecheck<<endl;
#endif
#endif
    double currentEnergy = FFM.computeEnergy(Bead::getDbData().coords.data(), Bead::getDbData().forces.data(), 0.0);
#ifdef DETAILEDOUTPUT_ENERGY
    CUDAcommon::handleerror(cudaDeviceSynchronize());
    double cuda_energy[1];
    CUDAcommon::handleerror(cudaMemcpy(cuda_energy, CUDAcommon::cudavars.gpu_energy,  sizeof(double),
                                       cudaMemcpyDeviceToHost));
    std::cout<<"Total Energy CE pN.nm CUDA "<<cuda_energy[0]<<" SERL "<<currentEnergy<<endl;
    std::cout<<endl;
#endif

    int iter = 0;
    while(!(cconvergencecheck[0])||!(sconvergencecheck)) {
        iter++;
        //TODO let each forcefield calculate energy IFF conv state = false. That will help
        // them avoid unnecessary iterations.
        //let each forcefield also add energies to two different energy variables.
        double energyLambda = FFM.computeEnergy(Bead::getDbData().coords.data(), Bead::getDbData().forces.data(), lambda);
#ifdef DETAILEDOUTPUT_ENERGY
        CUDAcommon::handleerror(cudaDeviceSynchronize());
        double cuda_energy[1];
        CUDAcommon::handleerror(cudaMemcpy(cuda_energy, CUDAcommon::cudavars.gpu_energy,  sizeof(double),
                                           cudaMemcpyDeviceToHost));
        std::cout<<"Total Energy EL pN.nm CUDA "<<cuda_energy[0]<<" SERL "
                ""<<energyLambda<<endl;
        std::cout<<endl;
#endif

#ifdef SERIAL
        //@{ Lambda phase 2
        if(!(sconvergencecheck)){
            double idealEnergyChange = -BACKTRACKSLOPE * lambda * allFDotFA();
            double energyChange = energyLambda - currentEnergy;
#ifdef DETAILEDOUTPUT_LAMBDA
            std::cout<<"BACKTRACKSLOPE "<<BACKTRACKSLOPE<<" lambda "<<lambda<<" allFDotFA"
                    " "<<allFDotFA()<<endl;
            std::cout<<"SL energyChange "<<energyChange<<" idealEnergyChange "
                    ""<<idealEnergyChange<<endl;
#endif
            //return if ok
            if(energyChange <= idealEnergyChange) {
                sconvergencecheck = true;}
            else
                //reduce lambda
                lambda *= LAMBDAREDUCE;

            if(lambda <= 0.0 || lambda <= LAMBDATOL) {
                sconvergencecheck = true;
                lambda = 0.0;

            }
#ifdef DETAILEDOUTPUT_LAMBDA
            std::cout<<"SL2 BACKTRACKSLOPE "<<BACKTRACKSLOPE<<" lambda "<<lambda<<" allFDotFA "
                                                                                <<allFDotFA()<<endl;
            std::cout<<"SL2 energyChange "<<energyChange<<" idealEnergyChange "
                    ""<<idealEnergyChange
                     <<" lambda "<<lambda<<" state "<<sconvergencecheck<<endl;
#endif
        }
        //@{ Lambda phase 2
#endif
    }
//    std::cout<<"lambda determined in "<<iter<< " iterations "<<endl;
//synchronize streams
    if(cconvergencecheck[0]||sconvergencecheck) {
#ifdef SERIAL
        delete [] cconvergencecheck;
#endif
        return lambda;
    }

}

double CGMethod::safeBacktrackingLineSearch(ForceFieldManager& FFM, double MAXDIST,
                                            double LAMBDAMAX, bool *gpu_safestate) {
    //reset safe mode
    _safeMode = false;
    sconvergencecheck = true;
    //calculate first lambda
    double lambda = LAMBDAMAX;
//    std::cout<<"safe 0"<<endl;
#ifdef SERIAL //SERIAL
    sconvergencecheck = false;
    cconvergencecheck = new bool[1];
    cconvergencecheck[0] = true;
#endif
//prepare for ping pong optimization
    double currentEnergy = FFM.computeEnergy(Bead::getDbData().coords.data(), Bead::getDbData().forces.data(), 0.0);
#ifdef DETAILEDOUTPUT_ENERGY
    CUDAcommon::handleerror(cudaDeviceSynchronize());
    double cuda_energy[1];
    CUDAcommon::handleerror(cudaMemcpy(cuda_energy, CUDAcommon::cudavars.gpu_energy,  sizeof(double),
                                       cudaMemcpyDeviceToHost));
    std::cout<<"Total Energy CE pN.nm CUDA "<<cuda_energy[0]<<" SERL "<<currentEnergy<<endl;
    std::cout<<endl;
#endif
    int iter =0;
    //safe backtracking loop
    while(!(cconvergencecheck[0])||!(sconvergencecheck)) {
        //new energy when moved by lambda
//        std::cout<<"safe z"<<endl;
        iter++;
        double energyLambda = FFM.computeEnergy(Bead::getDbData().coords.data(), Bead::getDbData().forces.data(), lambda);
#ifdef DETAILEDOUTPUT_ENERGY
        CUDAcommon::handleerror(cudaDeviceSynchronize());
        double cuda_energy[1];
        CUDAcommon::handleerror(cudaMemcpy(cuda_energy, CUDAcommon::cudavars.gpu_energy,  sizeof(double),
                                           cudaMemcpyDeviceToHost));
        std::cout<<"Total Energy EL pN.nm CUDA "<<cuda_energy[0]<<" SERL "
                ""<<energyLambda<<endl;
        std::cout<<endl;
#endif

#ifdef SERIAL
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
//            std::cout<<"safe energyChange "<<energyChange<<" lambda "<<lambda<<endl;
        }
#endif
    }
//    std::cout<<"lambda determined in "<<iter<< " iterations "<<endl;
    // FIXME: Make sure every branch returns
    if(cconvergencecheck[0]||sconvergencecheck) {
#ifdef SERIAL
        delete [] cconvergencecheck;
#endif
        return lambda;
    }
}
