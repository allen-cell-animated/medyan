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

__global__ void addvectorM(double *U, int *params, double *U_sum, double *U_tot){
  U_sum[0] = 0.0;
    double sum = 0.0;
  for(auto i=0;i<params[1];i++){
    if(U[i] == -1.0 && sum != -1.0){
        U_sum[0] = -1.0;
        U_tot[0] = -1.0;
        sum = -1.0;
      break;
    }
    else
      sum  += U[i];
  }
    U_sum[0] = sum;
    atomicAdd(&U_tot[0], sum);

}
__global__ void testifitworks(double *coord, double* c){
    const unsigned int thread_idx = (blockIdx.x * blockDim.x) + threadIdx.x;
    c[thread_idx] = 10.59;

}

__global__ void MotorGhostStretchingHarmonicenergy(double *coord, double *force, int *beadSet, double *kstr,
                                                   double *eql, double *pos1, double *pos2, int *params,
                                                   double *U_i) {

    extern __shared__ double s[];
    double *c1 = s;
    double *c2 = &c1[3 * blockDim.x];
    double *c3 = &c2[3 * blockDim.x];
    double *c4 = &c3[3 * blockDim.x];
    double v1[3], v2[3];
    double dist;
    int nint = params[1];
    int n = params[0];
    const unsigned int thread_idx = (blockIdx.x * blockDim.x) + threadIdx.x;

//    double *coord1, *coord2, *coord3, *coord4, dist, U_i;
//    double *v1 = new double[3];
//    double *v2 = new double[3];
    if(thread_idx<nint) {
        for(auto i=0;i<3;i++){
            U_i[thread_idx] =0.0;
            c1[3 * threadIdx.x + i] = coord[3 * beadSet[n * thread_idx] + i];
            c2[3 * threadIdx.x + i] = coord[3 * beadSet[n * thread_idx + 1] + i];
            c3[3 * threadIdx.x + i] = coord[3 * beadSet[n * thread_idx + 2] + i];
            c4[3 * threadIdx.x + i] = coord[3 * beadSet[n * thread_idx + 3] + i];
//            printf("%f %f %f %f \n", c1[3 * threadIdx.x + i],c2[3 * threadIdx.x + i],c3[3 * threadIdx.x + i],
//                   c4[3 * threadIdx.x + i]);
        }

    }
    __syncthreads();
//    printf("%i \n",nint);
    if(thread_idx<nint) {
//        printf("%f %f \n", pos1[thread_idx],pos2[thread_idx]);
        midPointCoordinate(v1, c1, c2, pos1[thread_idx], 3 * threadIdx.x);
//        printf("%f %f %f \n",v1[0],v1[1],v1[2]);
        midPointCoordinate(v2, c3, c4, pos2[thread_idx], 3 * threadIdx.x);
//        printf("%f %f %f \n",v2[0],v2[1],v2[2]);
        dist = twoPointDistance(v1, v2) -eql[thread_idx] ;
//        printf("%f \n",dist);
//        U_i[thread_idx] = dist;
        U_i[thread_idx] = 0.5 * kstr[thread_idx] * dist * dist;

        if (fabs(U_i[thread_idx]) == __longlong_as_double(0x7ff0000000000000) //infinity
            || U_i[thread_idx] != U_i[thread_idx] || U_i[thread_idx] < -1.0) {
            //TODO set culprit in host after return
//            MotorGhostInteractions::_motorCulprit = MotorGhost::getMotorGhosts()[i];
            U_i[thread_idx]=-1.0;
//            assert(0);
        }
//        printf("%f %f %f \n", dist, kstr[thread_idx], U_i[thread_idx]);
    }

//    __syncthreads();
}

__global__ void MotorGhostStretchingHarmonicenergyz(double *coord, double *f, int *beadSet, double *kstr,
                                                   double *eql, double *pos1, double *pos2, int *params,
                                                   double *U_i, double *z) {

    extern __shared__ double s[];
    double *c1 = s;
    double *c2 = &c1[3 * blockDim.x];
    double *c3 = &c2[3 * blockDim.x];
    double *c4 = &c3[3 * blockDim.x];
    double *f1 = &c4[3 * blockDim.x];
    double *f2 = &f1[3 * blockDim.x];
    double *f3 = &f2[3 * blockDim.x];
    double *f4 = &f3[3 * blockDim.x];

    double v1[3], v2[3];
    double dist;
    int nint = params[1];
    int n = params[0];
    const unsigned int thread_idx = (blockIdx.x * blockDim.x) + threadIdx.x;

//    printf("%d \n", nint);
    if(thread_idx<nint) {
        U_i[thread_idx] = 0.0;
        for(auto i=0;i<3;i++){
            c1[3 * threadIdx.x + i] = coord[3 * beadSet[n * thread_idx] + i];
            c2[3 * threadIdx.x + i] = coord[3 * beadSet[n * thread_idx + 1] + i];
            c3[3 * threadIdx.x + i] = coord[3 * beadSet[n * thread_idx + 2] + i];
            c4[3 * threadIdx.x + i] = coord[3 * beadSet[n * thread_idx + 3] + i];
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

//    if(thread_idx<nint) {
////        printf("%f \n", pos1[thread_idx]);
////        char ccc='a';
////        printf("%f %f %f %f %f %f %c \n",c1[3 * threadIdx.x + 0],c1[3 * threadIdx.x + 1],c1[3 * threadIdx.x + 2],f1[3
////                                                                                                                    *
////                                                                                                                            threadIdx.x + 0],f1[3 * threadIdx.x + 1],f1[3 * threadIdx.x + 2],ccc);
////        printf("%f %f %f %f %f %f \n",c2[3 * threadIdx.x + 0],c2[3 * threadIdx.x + 1],c2[3 * threadIdx.x + 2],f2[3 * threadIdx.x + 0],f2[3 * threadIdx.x + 1],f2[3 * threadIdx.x + 2]);
////        printf("%f \n",z[0]);
////
////        printf("%f %f %f %c\n",v1[0],v1[1],v1[2],ccc);
////
////        printf("%f %f %f %f %f %f %c \n",c3[3 * threadIdx.x + 0],c3[3 * threadIdx.x + 1],c3[3 * threadIdx.x + 2],f3[3
////                                                                                                                    *
////                                                                                                                    threadIdx.x + 0],f3[3 * threadIdx.x + 1],f3[3 * threadIdx.x + 2],ccc);
////        printf("%f %f %f %f %f %f \n",c4[3 * threadIdx.x + 0],c4[3 * threadIdx.x + 1],c4[3 * threadIdx.x + 2],f4[3 * threadIdx.x + 0],f4[3 * threadIdx.x + 1],f4[3 * threadIdx.x + 2]);
//
//        midPointCoordinateStretched(v1, c1, f1, c2, f2, pos1[thread_idx], z[0], 3 * threadIdx.x);
//        midPointCoordinateStretched(v2, c3, f3, c4, f4, pos2[thread_idx], z[0], 3 * threadIdx.x);
////        printf("%f %f %f \n",v2[0],v2[1],v2[2]);
//        __syncthreads();
////
//        dist = twoPointDistance(v1,  v2) - eql[thread_idx];
////        printf("%f \n", twoPointDistance(v1,  v2));
////        U_i[thread_idx] = 0.5 * kstr[thread_idx] * dist * dist;
//
//
////        U_i[thread_idx] =  kstr[thread_idx];
////        U_i[thread_idx] = c1[3 * threadIdx.x];
////        gU_i[thread_idx] = double(nint);
////        checkU[thread_idx] = beadSet[n * thread_idx];
//
//        if (fabs(U_i[thread_idx]) == __longlong_as_double(0x7ff0000000000000) //infinity
//            || U_i[thread_idx] != U_i[thread_idx] || U_i[thread_idx] < -1.0) {
//            //TODO set culprit in host after return
////            MotorGhostInteractions::_motorCulprit = MotorGhost::getMotorGhosts()[i];
//            U_i[thread_idx]=-1.0;
////            assert(0);
//        }
//
//    }
//
//    __syncthreads();

    if(thread_idx<nint) {
        midPointCoordinateStretched(v1, c1, f1, c2, f2, pos1[thread_idx], z[0], 3 * threadIdx.x);
        midPointCoordinateStretched(v2, c3, f3, c4, f4, pos2[thread_idx], z[0], 3 * threadIdx.x);

        dist = twoPointDistance(v1,  v2) - eql[thread_idx];
        U_i[thread_idx] = 0.5 * kstr[thread_idx] * dist * dist;

//        for(auto i=0;i<3;i++)
//            checkU[33 * threadIdx.x + i]=v1[i];
//        for(auto i=0;i<3;i++)
//            checkU[33 * threadIdx.x + 3 + i]=v2[i];
//        for(auto i=0;i<3;i++)
//            checkU[33 * threadIdx.x + 6 + i]=c1[3 * threadIdx.x + i];
//        for(auto i=0;i<3;i++)
//            checkU[33 * threadIdx.x + 9 + i]=c3[3 * threadIdx.x + i];
//        for(auto i=0;i<3;i++)
//            checkU[33 * threadIdx.x + 12 + i]=f1[3 * threadIdx.x + i];
//        for(auto i=0;i<3;i++)
//            checkU[33 * threadIdx.x + 15 + i]=f3[3 * threadIdx.x + i];
//        for(auto i=0;i<3;i++)
//            checkU[33 * threadIdx.x + 18 + i]=c2[3 * threadIdx.x + i];
//        for(auto i=0;i<3;i++)
//            checkU[33 * threadIdx.x + 21 + i]=c4[3 * threadIdx.x + i];
//        for(auto i=0;i<3;i++)
//            checkU[33 * threadIdx.x + 24 + i]=f2[3 * threadIdx.x + i];
//        for(auto i=0;i<3;i++)
//            checkU[33 * threadIdx.x + 27 + i]=f4[3 * threadIdx.x + i];
//        checkU[33 * threadIdx.x + 30]=pos1[thread_idx];
//        checkU[33 * threadIdx.x + 31]=pos2[thread_idx];
//        checkU[33 * threadIdx.x + 32]=z[0];

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
                                          double *kstr, double *eql, double *pos1, double *pos2, int *params
                                                  ){

    extern __shared__ double s[];
    double *c1 = s;
    double *c2 = &c1[3 * blockDim.x];
    double *c3 = &c2[3 * blockDim.x];
    double *c4 = &c3[3 * blockDim.x];
//    double *v1, *v2;
//    v1 = new double[3];
//    v2 = new double[3];
    double v1[3], v2[3];
    double dist, invL, f0, f1[3], f2[3], f3[3], f4[3];
    int nint = params[1];
    int n = params[0];
    const unsigned int thread_idx = (blockIdx.x * blockDim.x) + threadIdx.x;

//    double *coord1, *coord2, *coord3, *coord4, dist, U_i;
//    double *v1 = new double[3];
//    double *v2 = new double[3];
//    printf("%d \n", thread_idx );
    if(thread_idx<nint) {
        for(auto i=0;i<3;i++){
            c1[3 * threadIdx.x + i] = coord[3 * beadSet[n * thread_idx] + i];
            c2[3 * threadIdx.x + i] = coord[3 * beadSet[n * thread_idx + 1] + i];
            c3[3 * threadIdx.x + i] = coord[3 * beadSet[n * thread_idx + 2] + i];
            c4[3 * threadIdx.x + i] = coord[3 * beadSet[n * thread_idx + 3] + i];
        }

    }
    __syncthreads();

//    printf("%.14f %.14f %.14f %.14f %.14f %.14f %.14f %.14f %.14f %.14f %.14f %.14f \n",c1[3 * threadIdx.x] , c1[3 *
//                                                                                                                 threadIdx.x +1] , c1[3 *
//                                                                                                                    threadIdx
//            .x+2],  c2[3 * threadIdx.x] , c2[3 * threadIdx.x +1] , c2[3 * threadIdx.x +2], c3[3 * threadIdx.x] , c3[3 *
//                                                                                                                    threadIdx.x +1] , c3[3 * threadIdx.x +2], c4[3 * threadIdx.x] , c4[3 * threadIdx.x +1] , c4[3 * threadIdx.x +2]);

   // __syncthreads();

//    if(thread_idx<nint) {
//    printf("%i \n", 3*thread_idx);
//    }
//__syncthreads();
    if(thread_idx<nint) {
        v1[0] = c1[3 * threadIdx.x] * (1.0 - pos1[thread_idx]) + pos1[thread_idx] * c2[3 * threadIdx.x];
        v1[1] = c1[3 * threadIdx.x + 1] * (1.0 - pos1[thread_idx]) + pos1[thread_idx] * c2[3 * threadIdx.x + 1];
        v1[2] = c1[3 * threadIdx.x + 2] * (1.0 - pos1[thread_idx]) + pos1[thread_idx] * c2[3 * threadIdx.x + 2];

        midPointCoordinate(v1, c1, c2, pos1[thread_idx], 3 * threadIdx.x);
        // char ccc='a';
        //printf("%d %f %f %f \n", thread_idx,v1[0],v1[1],v1[2]);

        midPointCoordinate(v2, c3, c4, pos2[thread_idx], 3 * threadIdx.x);
        //printf("%d %f %f %f \n", thread_idx, v2[0],v2[1],v2[2]);
//        dist = twoPointDistance(v1, v2) - eql[thread_idx];
        dist = twoPointDistance(v1, v2);
//        printf("%f \n", dist);
        invL = 1 / dist;
        f0 = kstr[thread_idx] * (dist - eql[thread_idx]) * invL;
        //printf(" %d %f \n", thread_idx, f0);

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

//        U_i[36 * threadIdx.x] = c1[3 * threadIdx.x];
//        U_i[36 * threadIdx.x +1] = c1[3 * threadIdx.x +1];
//        U_i[36 * threadIdx.x +2] = c1[3 * threadIdx.x +2];
//        U_i[36 * threadIdx.x +3] = c2[3 * threadIdx.x];
//        U_i[36 * threadIdx.x +4] = c2[3 * threadIdx.x +1];
//        U_i[36 * threadIdx.x +5] = c2[3 * threadIdx.x +2];
//        U_i[36 * threadIdx.x +6] = c3[3 * threadIdx.x];
//        U_i[36 * threadIdx.x +7] = c3[3 * threadIdx.x +1];
//        U_i[36 * threadIdx.x +8] = c3[3 * threadIdx.x +2];
//        U_i[36 * threadIdx.x +9] = c4[3 * threadIdx.x];
//        U_i[36 * threadIdx.x +10] = c4[3 * threadIdx.x +1];
//        U_i[36 * threadIdx.x +11] = c4[3 * threadIdx.x +2];
//        U_i[36 * threadIdx.x +12] = f1[0];
//        U_i[36 * threadIdx.x +13] = f1[1];
//        U_i[36 * threadIdx.x +14] = f1[2];
//        U_i[36 * threadIdx.x +15] = f2[0];
//        U_i[36 * threadIdx.x +16] = f2[1];
//        U_i[36 * threadIdx.x +17] = f2[2];
//        U_i[36 * threadIdx.x +18] = f3[0];
//        U_i[36 * threadIdx.x +19] = f3[1];
//        U_i[36 * threadIdx.x +20] = f3[2];
//        U_i[36 * threadIdx.x +21] = f4[0];
//        U_i[36 * threadIdx.x +22] = f4[1];
//        U_i[36 * threadIdx.x +23] = f4[2];
//        for(auto i=0;i<3;i++)
//            U_i[36 * threadIdx.x +24 + i] = f[3 * beadSet[n * thread_idx] + i];
//        for(auto i=0;i<3;i++)
//            U_i[36 * threadIdx.x +27 + i] = f[3 * beadSet[n * thread_idx +1] + i];
//        for(auto i=0;i<3;i++)
//            U_i[36 * threadIdx.x +30 + i] = f[3 * beadSet[n * thread_idx +2] + i];
//        for(auto i=0;i<3;i++)
//            U_i[36 * threadIdx.x +33 + i] = f[3 * beadSet[n * thread_idx +3] + i];


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
