//
// Created by aravind on 1/25/18.
//

#ifndef CUDA_VEC_BRANCHINGSTRETCHINGHARMONICCUDA_H
#define CUDA_VEC_BRANCHINGSTRETCHINGHARMONICCUDA_H
#include "BranchingStretchingHarmonic.h"

#include "Bead.h"

#include "MathFunctions.h"
#include "SysParams.h"
#include <limits>
#include <assert.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include <cmath>
#include <math.h>

using namespace mathfunc;

#ifdef CUDAACCL

//#if !defined(__CUDA_ARCH__) || __CUDA_ARCH__ >= 600
//
//#else
//static __inline__ __device__ double atomicAdd(double *address, double val) {
//    unsigned long long int* address_as_ull = (unsigned long long int*)address;
//    unsigned long long int old = *address_as_ull, assumed;
//    if (val==0.0)
//      return __longlong_as_double(old);
//    do {
//      assumed = old;
//      old = atomicCAS(address_as_ull, assumed, __double_as_longlong(val +__longlong_as_double(assumed)));
//    } while (assumed != old);
//    return __longlong_as_double(old);
//  }
//
//
//#endif

//__global__ void addvectorBS(double *U, int *params, double *U_sum, double *U_tot){
//    U_sum[0] = 0.0;
//    double sum = 0.0;
//    for(auto i=0;i<params[1];i++){
//        if(U[i] == -1.0 && sum != -1.0){
//            U_sum[0] = -1.0;
//            U_tot[0] = -1.0;
//            sum = -1.0;
//            break;
//        }
//        else
//            sum  += U[i];
//    }
//    U_sum[0] = sum;
//    atomicAdd(&U_tot[0], sum);
//
//}

__global__ void BranchingStretchingHarmonicenergy(double *coord, double *force, int *beadSet, double *kstr,
                                                 double *eql, double *pos, int *params, double *U_i) {

    extern __shared__ double s[];
    double *c1 = s;
    double *c2 = &c1[3 * blockDim.x];
    double *c3 = &c2[3 * blockDim.x];
    double dist;
    double v1[3];
    int nint = params[1];
    int n = params[0];
    const unsigned int thread_idx = (blockIdx.x * blockDim.x) + threadIdx.x;

    if(thread_idx<nint) {
        for(auto i=0;i<3;i++){
            U_i[thread_idx] =0.0;
            c1[3 * threadIdx.x + i] = coord[3 * beadSet[n * thread_idx] + i];
            c2[3 * threadIdx.x + i] = coord[3 * beadSet[n * thread_idx + 1] + i];
            c3[3 * threadIdx.x + i] = coord[3 * beadSet[n * thread_idx + 2] + i];
        }

    }
    __syncthreads();
    if(thread_idx<nint) {
        midPointCoordinate(v1, c1, c2, pos[thread_idx], 3 * threadIdx.x);
        dist = twoPointDistancemixedID(v1, c3, 0, 3 * threadIdx.x) - eql[thread_idx];

        U_i[thread_idx] = 0.5 * kstr[thread_idx] * dist * dist;

        if (fabs(U_i[thread_idx]) == __longlong_as_double(0x7ff0000000000000) //infinity
            || U_i[thread_idx] != U_i[thread_idx] || U_i[thread_idx] < -1.0) {
            //TODO set culprit in host after return
//            FilamentInteractions::_motorCulprit = Filament::getFilaments()[i];
            U_i[thread_idx]=-1.0;
//            assert(0);
        }
//        printf("%f %f %f \n", dist, kstr[thread_idx], U_i[thread_idx]);
    }

//    __syncthreads();
}

__global__ void BranchingStretchingHarmonicenergyz(double *coord, double *f, int *beadSet, double *kstr,
                                                  double *eql, double *pos, int *params, double *U_i, double *z) {

    extern __shared__ double s[];
    double *c1 = s;
    double *c2 = &c1[3 * blockDim.x];
    double *c3 = &c2[3 * blockDim.x];
    double *f1 = &c3[3 * blockDim.x];
    double *f2 = &f1[3 * blockDim.x];
    double *f3 = &f2[3 * blockDim.x];

    double dist;
    double v1[3];
    double vzero[3]; vzero[0] = 0; vzero[1] = 0; vzero[2] = 0;
    int nint = params[1];
    int n = params[0];
    const unsigned int thread_idx = (blockIdx.x * blockDim.x) + threadIdx.x;

    if(thread_idx<nint) {
        U_i[thread_idx] = 0.0;
        for(auto i=0;i<3;i++){
            c1[3 * threadIdx.x + i] = coord[3 * beadSet[n * thread_idx] + i];
            c2[3 * threadIdx.x + i] = coord[3 * beadSet[n * thread_idx + 1] + i];
            c3[3 * threadIdx.x + i] = coord[3 * beadSet[n * thread_idx + 2] + i];
            f1[3 * threadIdx.x + i] = f[3 * beadSet[n * thread_idx] + i];
            f2[3 * threadIdx.x + i] = f[3 * beadSet[n * thread_idx + 1] + i];
            f3[3 * threadIdx.x + i] = f[3 * beadSet[n * thread_idx + 2] + i];
        }

    }
    __syncthreads();

    if(thread_idx<nint) {

        midPointCoordinateStretched(v1, c1, f1, c2, f2, pos[thread_idx], z[0], 3 * threadIdx.x);
        dist = twoPointDistanceStretchedmixedID(v1, vzero, c3, f3, z[0], 0, 3 * threadIdx.x) - eql[thread_idx];

        U_i[thread_idx] = 0.5 * kstr[thread_idx] * dist * dist;


        if (fabs(U_i[thread_idx]) == __longlong_as_double(0x7ff0000000000000) //infinity
            || U_i[thread_idx] != U_i[thread_idx] || U_i[thread_idx] < -1.0) {
            //TODO set culprit in host after return
//            FilamentInteractions::_motorCulprit = Filament::getFilaments()[i];
            U_i[thread_idx]=-1.0;
//            assert(0);
        }

    }

}


__global__ void BranchingStretchingHarmonicforces(double *coord, double *f, int *beadSet,
                                                 double *kstr, double *eql, double *pos, int *params){
    extern __shared__ double s[];
    double *c1 = s;
    double *c2 = &c1[3 * blockDim.x];
    double *c3 = &c2[3 * blockDim.x];
    double dist, invL, f0;
    double v1[3], f1[3], f2[3], f3[3];
    int nint = params[1];
    int n = params[0];
    const unsigned int thread_idx = (blockIdx.x * blockDim.x) + threadIdx.x;

    if(thread_idx<nint) {
        for(auto i=0;i<3;i++){
            c1[3 * threadIdx.x + i] = coord[3 * beadSet[n * thread_idx] + i];
            c2[3 * threadIdx.x + i] = coord[3 * beadSet[n * thread_idx + 1] + i];
            c3[3 * threadIdx.x + i] = coord[3 * beadSet[n * thread_idx + 2] + i];
        }

    }
    __syncthreads();

    if(thread_idx<nint) {
        midPointCoordinate(v1, c1, c2, pos[thread_idx], 3 * threadIdx.x);
        dist = twoPointDistancemixedID(v1, c3, 0, 3 * threadIdx.x);
        invL = 1 / dist;
        f0 = kstr[thread_idx] * ( dist - eql[thread_idx]) * invL;

        f1[0] =  -f0 * ( c3[3 * threadIdx.x] - v1[0] ) * (pos[thread_idx] - 1);
        f1[1] =  -f0 * ( c3[3 * threadIdx.x + 1] - v1[1] ) * (pos[thread_idx] - 1);
        f1[2] =  -f0 * ( c3[3 * threadIdx.x + 2] - v1[2] ) * (pos[thread_idx] - 1);

        // force i+1
        f2[0] =  f0 * ( c3[3 * threadIdx.x] - v1[0] ) * pos[thread_idx];
        f2[1] =  f0 * ( c3[3 * threadIdx.x + 1] - v1[1] ) * pos[thread_idx];
        f2[2] =  f0 * ( c3[3 * threadIdx.x + 2] - v1[2] ) * pos[thread_idx];

        //force on j
        f3[0] =  -f0 * ( c3[3 * threadIdx.x] - v1[0] );
        f3[1] =  -f0 * ( c3[3 * threadIdx.x + 1] - v1[1] );
        f3[2] =  -f0 * ( c3[3 * threadIdx.x + 2] - v1[2] );

        for (int i = 0; i < 3; i++) {
            atomicAdd(&f[3 * beadSet[n * thread_idx] + i], f1[i]);
            atomicAdd(&f[3 * beadSet[n * thread_idx + 1] + i], f2[i]);
            atomicAdd(&f[3 * beadSet[n * thread_idx + 2] + i], f3[i]);
        }
    }
}

#endif
#endif //CUDA_VEC_BRANCHINGSTRETCHINGHARMONICCUDA_H
