//
// Created by aravind on 9/29/17.
//

#ifndef CUDA_VEC_MOTORGHOSTSTRETCHINGHARMONICCUDA_H
#define CUDA_VEC_MOTORGHOSTSTRETCHINGHARMONICCUDA_H

#include "MotorGhostStretchingHarmonic.h"

#include "MotorGhostStretching.h"

#include "Bead.h"

#include "MathFunctions.h"
#include "SysParams.h"
#include <limits>
#include <assert.h>
#include <cuda.h>
#include <cuda_runtime.h>

using namespace mathfunc;

#ifdef CUDAACCL

#if !defined(__CUDA_ARCH__) || __CUDA_ARCH__ >= 600

#else
static __inline__ __device__ double atomicAdd(double *address, double val) {
    unsigned long long int* address_as_ull = (unsigned long long int*)address;
    unsigned long long int old = *address_as_ull, assumed;
    if (val==0.0)
      return __longlong_as_double(old);
    do {
      assumed = old;
      old = atomicCAS(address_as_ull, assumed, __double_as_longlong(val +__longlong_as_double(assumed)));
    } while (assumed != old);
    return __longlong_as_double(old);
  }


#endif

__global__ void testifitworks(double *coord, double* c){
    const unsigned int thread_idx = (blockIdx.x * blockDim.x) + threadIdx.x;
    c[thread_idx] = 10.59;

}

__global__ void MotorGhostStretchingHarmonicenergy(double *coord, double *force, int *beadSet, double *kstr,
                                                   double *eql, double *pos1, double *pos2, int *params,
                                                   double *U_i, double *gU_i, double *gc1, double *gc2, double *checkU) {

    extern __shared__ double s[];
    double *c1 = s;
    double *c2 = &c1[3 * blockDim.x];
    double *c3 = &c2[3 * blockDim.x];
    double *c4 = &c3[3 * blockDim.x];
    double *v1, *v2;
    v1 = new double[3];
    v2 = new double[3];
    double dist;
    int nint = params[1];
    int n = params[0];
    const unsigned int thread_idx = (blockIdx.x * blockDim.x) + threadIdx.x;

//    double *coord1, *coord2, *coord3, *coord4, dist, U_i;
//    double *v1 = new double[3];
//    double *v2 = new double[3];
    gU_i[thread_idx] = 9.23;
    if(thread_idx<nint) {
        for(auto i=0;i<3;i++){
            c1[3 * threadIdx.x + i] = coord[3 * beadSet[n * thread_idx] + i];
            c2[3 * threadIdx.x + i] = coord[3 * beadSet[n * thread_idx + 1] + i];
            c3[3 * threadIdx.x + i] = coord[3 * beadSet[n * thread_idx + 2] + i];
            c4[3 * threadIdx.x + i] = coord[3 * beadSet[n * thread_idx + 3] + i];
        }

    }
    __syncthreads();

    if(thread_idx<nint) {

        midPointCoordinate(v1, c1, c2, pos1[thread_idx], 3 * threadIdx.x);
        midPointCoordinate(v2, c3, c4, pos2[thread_idx], 3 * threadIdx.x);
        dist = twoPointDistance(v1, v2) -eql[thread_idx] ;
//        U_i[thread_idx] = dist;
        U_i[thread_idx] = 0.5 * kstr[thread_idx] * dist * dist;

        if (fabs(U_i[thread_idx]) == __longlong_as_double(0x7ff0000000000000) //infinity
            || U_i[thread_idx] != U_i[thread_idx] || U_i[thread_idx] < -1.0) {
            //TODO set culprit in host after return
//            MotorGhostInteractions::_motorCulprit = MotorGhost::getMotorGhosts()[i];
            U_i[thread_idx]=-1.0;
//            assert(0);
        }
    }

    __syncthreads();
}

__global__ void MotorGhostStretchingHarmonicenergyz(double *coord, double *f, int *beadSet, double *kstr,
                                                   double *eql, double *pos1, double *pos2, int *params,
                                                   double *U_i, double *gU_i, double *z, double *gc1, double *gc2, double *checkU) {

    extern __shared__ double s[];
    double *c1 = s;
    double *c2 = &c1[3 * blockDim.x];
    double *c3 = &c2[3 * blockDim.x];
    double *c4 = &c3[3 * blockDim.x];
    double *f1 = &c4[3 * blockDim.x];
    double *f2 = &f1[3 * blockDim.x];
    double *f3 = &f2[3 * blockDim.x];
    double *f4 = &f3[3 * blockDim.x];

    double *v1, *v2;
    v1 = new double[3];
    v2 = new double[3];
    double dist;
    int nint = params[1];
    int n = params[0];
    const unsigned int thread_idx = (blockIdx.x * blockDim.x) + threadIdx.x;

//    printf("%d \n", nint);
    if(thread_idx<nint) {
        for(auto i=0;i<3;i++){
            c1[3 * threadIdx.x + i] = coord[3 * beadSet[n * thread_idx] + i];
            c2[3 * threadIdx.x + i] = coord[3 * beadSet[n * thread_idx + 1] + i];
            c3[3 * threadIdx.x + i] = coord[3 * beadSet[n * thread_idx + 2] + i];
            c4[3 * threadIdx.x + i] = coord[3 * beadSet[n * thread_idx + 3] + i];
//            f1[3 * threadIdx.x + i] = f[3 * beadSet[n * thread_idx] + i];
//            f2[3 * threadIdx.x + i] = f[3 * beadSet[n * thread_idx + 1] + i];
//            f3[3 * threadIdx.x + i] = f[3 * beadSet[n * thread_idx + 2] + i];
//            f4[3 * threadIdx.x + i] = f[3 * beadSet[n * thread_idx + 3] + i];
        }

    }
    __syncthreads();
//
//        printf("%f %f %f %f %f %f %f %f %f %f %f %f \n", c1[3 * threadIdx.x ],c1[3 * threadIdx.x + 1],c1[3 *
//        threadIdx.x + 2],c2[3 * threadIdx.x],c2[3 * threadIdx.x + 1],c2[3 * threadIdx.x + 2],c3[3 * threadIdx.x],
//        c3[3 * threadIdx.x + 1],c3[3 * threadIdx.x + 2],c4[3 * threadIdx.x ],c4[3 * threadIdx.x + 1],c4[3 * threadIdx
//        .x + 2]);
//    __syncthreads();

    if(thread_idx<nint) {
        for(auto i=0;i<3;i++){
//            c1[3 * threadIdx.x + i] = coord[3 * beadSet[n * thread_idx] + i];
//            c2[3 * threadIdx.x + i] = coord[3 * beadSet[n * thread_idx + 1] + i];
//            c3[3 * threadIdx.x + i] = coord[3 * beadSet[n * thread_idx + 2] + i];
//            c4[3 * threadIdx.x + i] = coord[3 * beadSet[n * thread_idx + 3] + i];
            f1[3 * threadIdx.x + i] = f[3 * beadSet[n * thread_idx] + i];
            f2[3 * threadIdx.x + i] = f[3 * beadSet[n * thread_idx + 1] + i];
            f3[3 * threadIdx.x + i] = f[3 * beadSet[n * thread_idx + 2] + i];
            f4[3 * threadIdx.x + i] = f[3 * beadSet[n * thread_idx + 3] + i];
        }

    }
    __syncthreads();

//    printf("%f %f %f %f %f %f %f %f %f %f %f %f \n", f1[3 * threadIdx.x ],f1[3 * threadIdx.x + 1],f1[3 *
//                                                                                                     threadIdx.x +
//                   2],f2[3 * threadIdx.x],f2[3 * threadIdx.x + 1],f2[3 * threadIdx.x + 2],f3[3 * threadIdx.x],
//           f3[3 * threadIdx.x + 1],f3[3 * threadIdx.x + 2],f4[3 * threadIdx.x ],f4[3 * threadIdx.x + 1],f4[3 * threadIdx
//                    .x + 2]);
//    __syncthreads();

    if(thread_idx<nint) {
        printf("%f \n", pos1[thread_idx]);
//            dist = 5.0;
//        U_i[thread_idx] = 2.0;
        midPointCoordinateStretched(v1, c1, f1, c2, f2, pos1[thread_idx], z[0], 3 * threadIdx.x);
        midPointCoordinateStretched(v2, c3, f3, c4, f4, pos2[thread_idx], z[0], 3 * threadIdx.x);

        dist = twoPointDistance(v1,  v2) - eql[thread_idx];
        U_i[thread_idx] = 0.5 * kstr[thread_idx] * dist * dist;

//        printf("%f %f %f %f %f %f %f\n",v1[0],v1[1],v1[2],v2[0],v2[1],v2[2],dist);
//        U_i[thread_idx] =  kstr[thread_idx];
//        U_i[thread_idx] = c1[3 * threadIdx.x];
//        gU_i[thread_idx] = double(nint);
//        checkU[thread_idx] = beadSet[n * thread_idx];

        if (fabs(U_i[thread_idx]) == __longlong_as_double(0x7ff0000000000000) //infinity
            || U_i[thread_idx] != U_i[thread_idx] || U_i[thread_idx] < -1.0) {
            //TODO set culprit in host after return
//            MotorGhostInteractions::_motorCulprit = MotorGhost::getMotorGhosts()[i];
            U_i[thread_idx]=-1.0;
//            assert(0);
        }

    }

}


__global__ void MotorGhostStretchingHarmonicforces(double *coord, double *f, int *beadSet,
                                          double *kstr, double *eql, double *pos1, double *pos2, int *params,
                                                   double *U_i){

    extern __shared__ double s[];
    double *c1 = s;
    double *c2 = &c1[3 * blockDim.x];
    double *c3 = &c2[3 * blockDim.x];
    double *c4 = &c3[3 * blockDim.x];
    double *v1, *v2;
    v1 = new double[3];
    v2 = new double[3];
    double dist, invL, f0, f1[3], f2[3], f3[3], f4[3];
    int nint = params[1];
    int n = params[0];
    const unsigned int thread_idx = (blockIdx.x * blockDim.x) + threadIdx.x;

//    double *coord1, *coord2, *coord3, *coord4, dist, U_i;
//    double *v1 = new double[3];
//    double *v2 = new double[3];

    if(thread_idx<nint) {
        for(auto i=0;i<3;i++){
            c1[3 * threadIdx.x + i] = coord[3 * beadSet[n * thread_idx] + i];
            c2[3 * threadIdx.x + i] = coord[3 * beadSet[n * thread_idx + 1] + i];
            c3[3 * threadIdx.x + i] = coord[3 * beadSet[n * thread_idx + 2] + i];
            c4[3 * threadIdx.x + i] = coord[3 * beadSet[n * thread_idx + 3] + i];
        }

    }
    __syncthreads();


    if(thread_idx<nint) {
        midPointCoordinate(v1, c1, c2, pos1[thread_idx], 3 * threadIdx.x);
        midPointCoordinate(v2, c3, c4, pos2[thread_idx], 3 * threadIdx.x);
        dist = twoPointDistance(v1, v2) - eql[thread_idx];

        dist = twoPointDistance(v1, v2);
        invL = 1 / dist;
        f0 = kstr[thread_idx] * (dist - eql[thread_idx]) * invL;


        //force on i
        f1[0] = -f0 * (v1[0] - v2[0]) * (1 - pos1[thread_idx]);
        f1[1] = -f0 * (v1[1] - v2[1]) * (1 - pos1[thread_idx]);
        f1[2] = -f0 * (v1[2] - v2[2]) * (1 - pos1[thread_idx]);

        // force i+1
        f2[0] = -f0 * (v1[0] - v2[0]) * (pos1[thread_idx]);
        f2[1] = -f0 * (v1[1] - v2[1]) * (pos1[thread_idx]);
        f2[2] = -f0 * (v1[2] - v2[2]) * (pos1[thread_idx]);

        //force on j
        f3[0] = f0 * (v1[0] - v2[0]) * (1 - pos2[thread_idx]);
        f3[1] = f0 * (v1[1] - v2[1]) * (1 - pos2[thread_idx]);
        f3[2] = f0 * (v1[2] - v2[2]) * (1 - pos2[thread_idx]);

        // force j+1
        f4[0] = f0 * (v1[0] - v2[0]) * (pos2[thread_idx]);
        f4[1] = f0 * (v1[1] - v2[1]) * (pos2[thread_idx]);
        f4[2] = f0 * (v1[2] - v2[2]) * (pos2[thread_idx]);

//        for (int i = 0; i < 3; i++) {
//            f[3 * beadSet[n * thread_idx] + i] = 1.2345;
//        }

        for (int i = 0; i < 3; i++) {
            atomicAdd(&f[3 * beadSet[n * thread_idx] + i], f1[i]);
            atomicAdd(&f[3 * beadSet[n * thread_idx + 1] + i], f2[i]);
            atomicAdd(&f[3 * beadSet[n * thread_idx + 2] + i], f3[i]);
            atomicAdd(&f[3 * beadSet[n * thread_idx + 3] + i], f4[i]);
        }
    }
    }

#endif
#endif //CUDA_VEC_MOTORGHOSTSTRETCHINGHARMONICCUDA_H
