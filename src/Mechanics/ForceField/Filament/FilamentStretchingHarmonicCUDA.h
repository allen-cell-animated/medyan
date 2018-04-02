//
// Created by aravind on 1/22/18.
//

#ifndef CUDA_VEC_FILAMENTSTRETCHINGHARMONICCUDA_H
#define CUDA_VEC_FILAMENTSTRETCHINGHARMONICCUDA_H
#include "FilamentStretchingHarmonic.h"

#include "FilamentStretching.h"

#include "Bead.h"

#include "MathFunctions.h"
#include "SysParams.h"
#include <limits>
#include <assert.h>
#include <cuda.h>
#include <cuda_runtime.h>

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

//__global__ void addvectorFS(double *U, int *params, double *U_sum, double *U_tot){
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

__global__ void FilamentStretchingHarmonicenergy(double *coord, double *force, int *beadSet, double *kstr,
                                                 double *eql, int *params, double *U_i, double *z, int *culpritID,
                                                 char* culpritFF, char* culpritinteraction, char* FF, char*
                                                 interaction) {
    if(z[0] == 0.0) {
        extern __shared__ double s[];
        double *c1 = s;
        double *c2 = &c1[3 * blockDim.x];
//    double v1[3], v2[3];
        double dist;
        int nint = params[1];
        int n = params[0];
        const unsigned int thread_idx = (blockIdx.x * blockDim.x) + threadIdx.x;

        if (thread_idx < nint) {
            for (auto i = 0; i < 3; i++) {
                U_i[thread_idx] = 0.0;
                c1[3 * threadIdx.x + i] = coord[3 * beadSet[n * thread_idx] + i];
                c2[3 * threadIdx.x + i] = coord[3 * beadSet[n * thread_idx + 1] + i];
            }

//        }
//        __syncthreads();
//        if (thread_idx < nint) {
            dist = twoPointDistance(c1, c2, 3 * threadIdx.x) - eql[thread_idx];
//        if(thread_idx ==1)
            U_i[thread_idx] = 0.5 * kstr[thread_idx] * dist * dist;
//        U_i[thread_idx] = __longlong_as_double(0x7ff0000000000000);
            if (fabs(U_i[thread_idx]) == __longlong_as_double(0x7ff0000000000000) //infinity
                || U_i[thread_idx] != U_i[thread_idx] || U_i[thread_idx] < -1.0) {

                U_i[thread_idx] = -1.0;
                culpritID[0] = thread_idx;
                culpritID[1] = -1;
                int j = 0;
                while (FF[j] != 0) {
                    culpritFF[j] = FF[j];
                    j++;
                }
                j = 0;
                while (interaction[j] != 0) {
                    culpritinteraction[j] = interaction[j];
                    j++;
                }
                assert(0);
                __syncthreads();
            }
        }
    }
//    __syncthreads();
}

__global__ void FilamentStretchingHarmonicenergyz(double *coord, double *f, int *beadSet, double *kstr,
                                                  double *eql, int *params, double *U_i, double *z, int *culpritID,
                                                  char* culpritFF, char* culpritinteraction, char* FF, char*
                                                  interaction) {
    if(z[0] != 0.0) {
        extern __shared__ double s[];
        double *c1 = s;
        double *c2 = &c1[3 * blockDim.x];
        double *f1 = &c2[3 * blockDim.x];
        double *f2 = &f1[3 * blockDim.x];

        //    double v1[3], v2[3];
        double dist;
        int nint = params[1];
        int n = params[0];
        const unsigned int thread_idx = (blockIdx.x * blockDim.x) + threadIdx.x;

        if (thread_idx < nint) {
            U_i[thread_idx] = 0.0;
            for (auto i = 0; i < 3; i++) {
                c1[3 * threadIdx.x + i] = coord[3 * beadSet[n * thread_idx] + i];
                c2[3 * threadIdx.x + i] = coord[3 * beadSet[n * thread_idx + 1] + i];
                f1[3 * threadIdx.x + i] = f[3 * beadSet[n * thread_idx] + i];
                f2[3 * threadIdx.x + i] = f[3 * beadSet[n * thread_idx + 1] + i];
            }
//     }
//        __syncthreads();
//
//        if (thread_idx < nint) {
           dist = twoPointDistanceStretched(c1, f1, c2, f2, z[0], 3 * threadIdx.x) - eql[thread_idx];
           U_i[thread_idx] = 0.5 * kstr[thread_idx] * dist * dist;
//        printf("C %d %f %f\n", thread_idx, dist, U_i[thread_idx]);

            if (fabs(U_i[thread_idx]) == __longlong_as_double(0x7ff0000000000000) //infinity
                || U_i[thread_idx] != U_i[thread_idx] || U_i[thread_idx] < -1.0) {
                U_i[thread_idx] = -1.0;
                culpritID[0] = thread_idx;
                culpritID[1] = -1;
                 int j = 0;
              while (FF[j] != 0) {
                    culpritFF[j] = FF[j];
                 j++;
             }
                j = 0;
                while (interaction[j] != 0) {
                    culpritinteraction[j] = interaction[j];
                    j++;
                }
                assert(0);
                __syncthreads();
            }

        }
    }
}


__global__ void FilamentStretchingHarmonicforces(double *coord, double *f, int *beadSet,
                                                   double *kstr, double *eql, int *params){
    extern __shared__ double s[];
    double *c1 = s;
    double *c2 = &c1[3 * blockDim.x];
    double dist, invL, f0, f1[3], f2[3];
    int nint = params[1];
    int n = params[0];
    const unsigned int thread_idx = (blockIdx.x * blockDim.x) + threadIdx.x;

    if(thread_idx<nint) {
        for(auto i=0;i<3;i++){
            c1[3 * threadIdx.x + i] = coord[3 * beadSet[n * thread_idx] + i];
            c2[3 * threadIdx.x + i] = coord[3 * beadSet[n * thread_idx + 1] + i];
//            printf("before %d %f %f %f %f %f %f\n", thread_idx, f[3 * beadSet[n * thread_idx] + 0],f[3 * beadSet[n * thread_idx] + 1],f[3 * beadSet[n *
//                                                                                                              thread_idx] + 2], f[3 * beadSet[n * thread_idx + 1]], f[3 * beadSet[n * thread_idx + 1]+1],
//                   f[3 *
//                     beadSet[n * thread_idx + 1]+2]);
        }
//    }
//    __syncthreads();
//
//    if(thread_idx<nint) {
//        v1[0] = c1[3 * threadIdx.x] * (1.0 - pos1[thread_idx]) + pos1[thread_idx] * c2[3 * threadIdx.x];
//        v1[1] = c1[3 * threadIdx.x + 1] * (1.0 - pos1[thread_idx]) + pos1[thread_idx] * c2[3 * threadIdx.x + 1];
//        v1[2] = c1[3 * threadIdx.x + 2] * (1.0 - pos1[thread_idx]) + pos1[thread_idx] * c2[3 * threadIdx.x + 2];
//
//        midPointCoordinate(v1, c1, c2, pos1[thread_idx], 3 * threadIdx.x);
//
//        midPointCoordinate(v2, c3, c4, pos2[thread_idx], 3 * threadIdx.x);
        dist = twoPointDistance(c1, c2, 3 * threadIdx.x);
        invL = 1 / dist;
        f0 = kstr[thread_idx] * (dist - eql[thread_idx]) * invL;

        //force on i
        f2[0] = f0 * (c1[3 * threadIdx.x] - c2[3 * threadIdx.x]);
        f2[1] = f0 * (c1[3 * threadIdx.x + 1] - c2[3 * threadIdx.x + 1]);
        f2[2] = f0 * (c1[3 * threadIdx.x + 2] - c2[3 * threadIdx.x + 2]);

        // force i+1
        f1[0] = f0 * (-c1[3 * threadIdx.x] + c2[3 * threadIdx.x]);
        f1[1] = f0 * (-c1[3 * threadIdx.x + 1] + c2[3 * threadIdx.x + 1]);
        f1[2] = f0 * (-c1[3 * threadIdx.x + 2] + c2[3 * threadIdx.x + 2]);
//        printf("Filcalc %d %d %d %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n", beadSet[n * thread_idx], beadSet[n * thread_idx +
//                       1],
//               thread_idx,dist, eql[thread_idx], kstr[thread_idx], f0,
//               f1[0], f1[1], f1[2], f2[0],
//               f2[1],
//               f2[2], c1[3 * threadIdx.x],c1[3 * threadIdx.x +1],c1[3 * threadIdx.x +2], c2[3 * threadIdx.x],c2[3 *
//               threadIdx.x +1],c2[3 * threadIdx.x +2]);

//        printf("%d %f %f %f %f %f %f %f %f\n", thread_idx, kstr[thread_idx], f0, c1[3 * threadIdx.x],c1[3 * threadIdx
//                                                                                                                    .x +1],c1[3 * threadIdx.x +2],c2[3 * threadIdx.x],c2[3 * threadIdx
//                .x +1],c2[3 * threadIdx.x +2]);
//        printf("%d %f %f %f %f %f %d %f %f %f %f %f %f\n", thread_idx, kstr[thread_idx], dist, eql[thread_idx], invL,
//               f0, f1[0],f1[1],f1[2],f2[0],f2[1],f2[2]);
//        printf("%d %f %f %f %f %f %f\n", thread_idx, f1[0],f1[1],f1[2],f2[0],f2[1],f2[2]);
        for (int i = 0; i < 3; i++) {
            if (fabs(f1[i]) == __longlong_as_double(0x7ff0000000000000) //infinity
                || f1[i] != f1[i]) {
                printf("Fil. Stret. Force became infinite %f %f %f\n",f1[0], f1[1], f1[2]);
                assert(0);
            }
            if (fabs(f2[i]) == __longlong_as_double(0x7ff0000000000000) //infinity
                || f2[i] != f2[i]) {
                printf("Fil. Stret. became infinite %f %f %f\n",f2[0], f2[1], f2[2]);
                assert(0);
            }
            atomicAdd(&f[3 * beadSet[n * thread_idx] + i], f1[i]);
            atomicAdd(&f[3 * beadSet[n * thread_idx + 1] + i], f2[i]);
        }

    }
}

#endif
#endif //CUDA_VEC_FILAMENTSTRETCHINGHARMONICCUDA_H
