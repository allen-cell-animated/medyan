//
// Created by aravind on 9/26/17.
//

#ifndef CUDA_VEC_CYLINDEREXCLVOLREPULSIONCUDA_H
#define CUDA_VEC_CYLINDEREXCLVOLREPULSIONCUDA_H
#include "CylinderExclVolRepulsion.h"

#include "CylinderExclVolume.h"

#include "Bead.h"

#include "MathFunctions.h"
#include "SysParams.h"
#include <limits>

#include <cuda.h>
#include <cuda_runtime.h>
#include <assert.h>
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

__global__
void saxpy(int n, float a, float *x, float *y)
{
    int i = blockIdx.x*blockDim.x + threadIdx.x;
    if (i < n) y[i] = a*x[i] + y[i];
}
__global__ void CUDAExclVolRepulsionenergy(double *coord, double *force, int *beadSet, double *krep,
                                           int *params, double *U_i, double *z, int *culpritID,
                                           char* culpritFF, char* culpritinteraction, char* FField, char*
                                           interaction) {
//memory needed: 34*THREADSPERBLOCK*sizeof(double)+2*THREADSPERBLOCK*sizeof(int);
    if(z[0] == 0.0) {
        extern __shared__ double s[];
//    double *coord_image=s;
//    double *c1 = &coord_image[3 * blockDim.x];
        double *c1 = s;
        double *c2 = &c1[3 * blockDim.x];
        double *c3 = &c2[3 * blockDim.x];
        double *c4 = &c3[3 * blockDim.x];
        double d, invDSquare;
        double a, b, c, e, F, AA, BB, CC, DD, EE, FF, GG, HH, JJ;
        double ATG1, ATG2, ATG3, ATG4;
        int nint = params[1];
        int n = params[0];
        const unsigned int thread_idx = (blockIdx.x * blockDim.x) + threadIdx.x;

        if (thread_idx < nint) {
            U_i[thread_idx] = 0.0;
            for (auto i = 0; i < 3; i++) {
                c1[3 * threadIdx.x + i] = coord[3 * beadSet[n * thread_idx] + i];
                c2[3 * threadIdx.x + i] = coord[3 * beadSet[n * thread_idx + 1] + i];
                c3[3 * threadIdx.x + i] = coord[3 * beadSet[n * thread_idx + 2] + i];
                c4[3 * threadIdx.x + i] = coord[3 * beadSet[n * thread_idx + 3] + i];
            }

//        }
//        __syncthreads();
//
//        if (thread_idx < nint) {
            //check if parallel
            if (areParallel(c1, c2, c3, c4, 3 * threadIdx.x)) {

                d = twoPointDistance(c1, c3, 3 * threadIdx.x);
                invDSquare = 1 / (d * d);
                U_i[thread_idx] = krep[thread_idx] * invDSquare * invDSquare;

                if (U_i[thread_idx] == __longlong_as_double(0x7ff0000000000000) //infinity
                    || U_i[thread_idx] != U_i[thread_idx] || U_i[thread_idx] < -1.0) {
                    U_i[thread_idx] = -1.0;
                    culpritID[0] = beadSet[n * thread_idx];
                    culpritID[1] = beadSet[n * thread_idx + 1];
                    culpritID[2] = beadSet[n * thread_idx + 2];
                    culpritID[3] = beadSet[n * thread_idx + 3];
                    int j = 0;
                    while (FField[j] != 0) {
                        culpritFF[j] = FField[j];
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

            } else {

                //check if in same plane
                if (areInPlane(c1, c2, c3, c4, 3 * threadIdx.x)) {

                    //slightly move point
                    movePointOutOfPlane(c1, c2, c3, c4, 2, 0.01, 3 * threadIdx.x);
                }

                a = scalarProduct(c1, c2, c1, c2, 3 * threadIdx.x);
                b = scalarProduct(c3, c4, c3, c4, 3 * threadIdx.x);
                c = scalarProduct(c3, c1, c3, c1, 3 * threadIdx.x);
                d = scalarProduct(c1, c2, c3, c4, 3 * threadIdx.x);
                e = scalarProduct(c1, c2, c3, c1, 3 * threadIdx.x);
                F = scalarProduct(c3, c4, c3, c1, 3 * threadIdx.x);

                AA = sqrt(a * c - e * e);
                BB = sqrt(b * c - F * F);

                CC = d * e - a * F;
                DD = b * e - d * F;

                EE = sqrt(a * (b + c - 2 * F) - (d - e) * (d - e));
                FF = sqrt(b * (a + c + 2 * e) - (d + F) * (d + F));

                GG = d * d - a * b - CC;
                HH = CC + GG - DD;
                JJ = c * (GG + CC) + e * DD - F * CC;


                ATG1 = atan((a + e) / AA) - atan(e / AA);
                ATG2 = atan((a + e - d) / EE) - atan((e - d) / EE);
                ATG3 = atan((F) / BB) - atan((F - b) / BB);
                ATG4 = atan((d + F) / FF) - atan((d + F - b) / FF);

                U_i[thread_idx] = 0.5 * krep[thread_idx] / JJ *
                                  (CC / AA * ATG1 + GG / EE * ATG2 + DD / BB * ATG3 + HH / FF * ATG4);
//        U_i[thread_idx] =a;
//        gU_i[thread_idx] = U_i[thread_idx];
//        for(auto j=0; j<3; j++) {
//            gc1[3 * threadIdx.x + j] = c1[3 * threadIdx.x + j];
//            gc2[3 * threadIdx.x + j] = c2[3 * threadIdx.x + j];
//        }

                if (U_i[thread_idx] == __longlong_as_double(0x7ff0000000000000)
                    || U_i[thread_idx] != U_i[thread_idx] || U_i[thread_idx] < -1.0) {
                    U_i[thread_idx] = -1.0;
                    culpritID[0] = beadSet[n * thread_idx];
                    culpritID[1] = beadSet[n * thread_idx + 1];
                    culpritID[2] = beadSet[n * thread_idx + 2];
                    culpritID[3] = beadSet[n * thread_idx + 3];
                    int j = 0;
                    while (FField[j] != 0) {
                        culpritFF[j] = FField[j];
                        j++;
                    }
                    j = 0;
                    while (interaction[j] != 0) {
                        culpritinteraction[j] = interaction[j];
                        j++;
                    }
                    printf("Coordiantes \n %f %f %f \n %f %f %f \n %f %f %f \n %f %f %f \n", c1[3 * threadIdx.x ],
                           c1[3 * threadIdx
                    .x + 1], c1[3 * threadIdx.x + 2], c2[3 * threadIdx.x ], c2[3 * threadIdx
                            .x + 1], c2[3 * threadIdx.x + 2], c3[3 * threadIdx.x ], c3[3 * threadIdx
                            .x + 1], c3[3 * threadIdx.x + 2], c4[3 * threadIdx.x ], c4[3 * threadIdx
                            .x + 1], c4[3 * threadIdx.x + 2]);
                    assert(0);
                    __syncthreads();
                }
            }
        }
//    else {
//        U_i[thread_idx] = 0.0;
////        gU_i[thread_idx] = U_i[thread_idx];
////        for(auto j=0; j<3; j++) {
////            gc1[3 * threadIdx.x + j] = 0.0;
////            gc2[3 * threadIdx.x + j] = 0.0;
////        }
//    }
//    printf("E %.16f \n", U_i[thread_idx] );
//    checkU[threadIdx.x]=0.0;
//    __syncthreads();

//    for(int offset = blockDim.x / 2;
//        offset > blockDim.x / 4; offset >>= 1)
//    {
//        if(threadIdx.x < offset)
//        {
////            if(U_i[threadIdx.x + offset]==-1.0||U_i[threadIdx.x]==-1.0)
////                U_i[threadIdx.x]=-1;
////            //add a partial sum upstream to our own
////            else
//            U_i[threadIdx.x] += U_i[threadIdx.x + offset];
//            checkU[threadIdx.x]=U_i[threadIdx.x];
////            gU_i[threadIdx.x]=U_i[threadIdx.x];
////                U_i[threadIdx.x]=offset;
//        }
////        if(threadIdx.x==0) {
////            int i=0;
////            int frac=offset;
////            while(frac!=1)
////            {frac=frac/2;i=i+1;}
////            checkU[threadIdx.x + i] = U_i[threadIdx.x];
////        }
//
//        // wait until all threads in the block have updated their partial sums
//        __syncthreads();
//
//    }
//
//    __syncthreads();
//    // thread 0 writes the final result
//    if(threadIdx.x == 0)
//    {
//        U[blockIdx.x] = U_i[0];
//    }
//    __syncthreads();
    }
}


__global__ void CUDAExclVolRepulsionenergyz(double *coord, double *f, int *beadSet,
                                            double *krep, int *params, double *U_i, double *z, int *culpritID,
                                            char* culpritFF, char* culpritinteraction, char* FField, char*
                                            interaction) {
    if(z[0] != 0.0) {

        extern __shared__ double s[];

        double *c1 = s;
        double *c2 = &c1[3 * blockDim.x];
        double *c3 = &c2[3 * blockDim.x];
        double *c4 = &c3[3 * blockDim.x];


        double d, invDSquare;
        double a, b, c, e, F, AA, BB, CC, DD, EE, FF, GG, HH, JJ;
        double ATG1, ATG2, ATG3, ATG4;
        int nint = params[1];
        int n = params[0];
        const unsigned int thread_idx = (blockIdx.x * blockDim.x) + threadIdx.x;

        if (thread_idx < nint) {
            U_i[thread_idx] = 0.0;
            for (auto i = 0; i < 3; i++) {
                c1[3 * threadIdx.x + i] = coord[3 * beadSet[n * thread_idx] + i]
                                          + z[0] * f[3 * beadSet[n * thread_idx] + i];
                c2[3 * threadIdx.x + i] = coord[3 * beadSet[n * thread_idx + 1] + i]
                                          + z[0] * f[3 * beadSet[n * thread_idx + 1] + i];
                c3[3 * threadIdx.x + i] = coord[3 * beadSet[n * thread_idx + 2] + i]
                                          + z[0] * f[3 * beadSet[n * thread_idx + 2] + i];
                c4[3 * threadIdx.x + i] = coord[3 * beadSet[n * thread_idx + 3] + i]
                                          + z[0] * f[3 * beadSet[n * thread_idx + 3] + i];
            }

//        }
//        __syncthreads();
//
//        if (thread_idx < nint) {
            //check if parallel
            if (areParallel(c1, c2, c3, c4, 3 * threadIdx.x)) {

                d = twoPointDistance(c1, c3, 3 * threadIdx.x);
                invDSquare = 1 / (d * d);
                U_i[thread_idx] = krep[thread_idx] * invDSquare * invDSquare;

                if (U_i[thread_idx] == __longlong_as_double(0x7ff0000000000000) //infinity
                    || U_i[thread_idx] != U_i[thread_idx] || U_i[thread_idx] < -1.0) {
                    U_i[thread_idx] = -1.0;
                    culpritID[0] = beadSet[n * thread_idx];
                    culpritID[1] = beadSet[n * thread_idx + 1];
                    culpritID[2] = beadSet[n * thread_idx + 2];
                    culpritID[3] = beadSet[n * thread_idx + 3];
                    int j = 0;
                    while (FField[j] != 0) {
                        culpritFF[j] = FField[j];
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
            } else {
//            auto jj=3;
                //check if in same plane
                if (areInPlane(c1, c2, c3, c4, 3 * threadIdx.x)) {
//            jj=2;
                    //slightly move point
                    movePointOutOfPlane(c1, c2, c3, c4, 2, 0.01, 3 * threadIdx.x);
//            printf("%i %i %.16f %.16f %.16f %.16f %.16f %.16f %.16f %.16f %.16f %.16f %.16f %.16f \n",thread_idx, jj,
//                           c1[3 * threadIdx.x], c1[3 *
//                                                                                                                                       threadIdx.x + 1], c1[3 * threadIdx.x + 2]
//                    , c2[3 * threadIdx.x], c2[3 * threadIdx.x + 1], c2[3 * threadIdx.x + 2], c3[3 * threadIdx.x], c3[3 * threadIdx.x + 1],
//                   c3[3 * threadIdx.x + 2], c4[3 * threadIdx.x], c4[3 * threadIdx.x + 1], c4[3 * threadIdx.x + 2]);
                }
                // else
//            printf("%i %i %.16f %.16f %.16f %.16f %.16f %.16f %.16f %.16f %.16f %.16f %.16f %.16f\n"
//                           ,thread_idx, jj, c1[3 * threadIdx.x], c1[3 *
//                                                                                                                                       threadIdx.x + 1], c1[3 * threadIdx.x + 2]
//                    , c2[3 * threadIdx.x], c2[3 * threadIdx.x + 1], c2[3 * threadIdx.x + 2], c3[3 * threadIdx.x], c3[3 * threadIdx.x + 1],
//                   c3[3 * threadIdx.x + 2], c4[3 * threadIdx.x], c4[3 * threadIdx.x + 1], c4[3 * threadIdx.x + 2]);
                a = scalarProduct(c1, c2, c1, c2, 3 * threadIdx.x);
                b = scalarProduct(c3, c4, c3, c4, 3 * threadIdx.x);
                c = scalarProduct(c3, c1, c3, c1, 3 * threadIdx.x);
                d = scalarProduct(c1, c2, c3, c4, 3 * threadIdx.x);
                e = scalarProduct(c1, c2, c3, c1, 3 * threadIdx.x);
                F = scalarProduct(c3, c4, c3, c1, 3 * threadIdx.x);

                AA = sqrt(a * c - e * e);
                BB = sqrt(b * c - F * F);

                CC = d * e - a * F;
                DD = b * e - d * F;

                EE = sqrt(a * (b + c - 2 * F) - (d - e) * (d - e));
                FF = sqrt(b * (a + c + 2 * e) - (d + F) * (d + F));

                GG = d * d - a * b - CC;
                HH = CC + GG - DD;
                JJ = c * (GG + CC) + e * DD - F * CC;


                ATG1 = atan((a + e) / AA) - atan(e / AA);
                ATG2 = atan((a + e - d) / EE) - atan((e - d) / EE);
                ATG3 = atan((F) / BB) - atan((F - b) / BB);
                ATG4 = atan((d + F) / FF) - atan((d + F - b) / FF);

                U_i[thread_idx] = 0.5 * krep[thread_idx] / JJ *
                                  (CC / AA * ATG1 + GG / EE * ATG2 + DD / BB * ATG3 + HH / FF * ATG4);
//        U_i[thread_idx] =a;
//        for(auto j=0; j<3; j++) {
//            gc1[3 * threadIdx.x + j] = f[3 * beadSet[n * thread_idx] + j];
//            gc2[3 * threadIdx.x + j] = f[3 * beadSet[n * thread_idx + 1] + j];
//        }

                if (U_i[thread_idx] == __longlong_as_double(0x7ff0000000000000)
                    || U_i[thread_idx] != U_i[thread_idx] || U_i[thread_idx] < -1.0) {
                    U_i[thread_idx] = -1.0;
                    culpritID[0] = beadSet[n * thread_idx];
                    culpritID[1] = beadSet[n * thread_idx + 1];
                    culpritID[2] = beadSet[n * thread_idx + 2];
                    culpritID[3] = beadSet[n * thread_idx + 3];
                    int j = 0;
                    while (FField[j] != 0) {
                        culpritFF[j] = FField[j];
                        j++;
                    }
                    j = 0;
                    while (interaction[j] != 0) {
                        culpritinteraction[j] = interaction[j];
                        j++;
                    }
                    printf("Coordiantes \n %f %f %f \n %f %f %f \n %f %f %f \n %f %f %f \n", c1[3 * threadIdx.x ],
                           c1[3 * threadIdx
                                   .x + 1], c1[3 * threadIdx.x + 2], c2[3 * threadIdx.x ], c2[3 * threadIdx
                                    .x + 1], c2[3 * threadIdx.x + 2], c3[3 * threadIdx.x ], c3[3 * threadIdx
                                    .x + 1], c3[3 * threadIdx.x + 2], c4[3 * threadIdx.x ], c4[3 * threadIdx
                                    .x + 1], c4[3 * threadIdx.x + 2]);
                    printf("Forces \n %f %f %f \n %f %f %f \n %f %f %f \n %f %f %f \n", c1[3 * threadIdx.x ],
                           f[3 * beadSet[n * thread_idx]], f[3 * beadSet[n * thread_idx] + 1], f[3 * beadSet[n *
                           thread_idx] + 2], f[3 * beadSet[n * thread_idx +1]], f[3 * beadSet[n * thread_idx +1] + 1],
                           f[3 * beadSet[n * thread_idx +1] + 2], f[3 * beadSet[n * thread_idx +2]], f[3 * beadSet[n *
                           thread_idx +2] + 1], f[3 * beadSet[n * thread_idx +2] + 2], f[3 * beadSet[n * thread_idx +3]],
                           f[3 * beadSet[n * thread_idx+3] + 1], f[3 * beadSet[n * thread_idx+3] + 2]);
                    assert(0);
                    __syncthreads();
                }
            }
        }
//    else {
//        U_i[thread_idx] = 0.0;
////        for(auto j=0; j<3; j++) {
////            gc1[3 * threadIdx.x + j] = 0.0;
////            gc2[3 * threadIdx.x + j] = 0.0;
////        }
//    }
        //printf("%.16f \n", U_i[thread_idx]);

//    __syncthreads();

//    for(int offset = blockDim.x / 2;
//        offset > 0; offset >>= 1)
//    {
//        if(threadIdx.x < offset)
//        {
////            if(U_i[threadIdx.x + offset]==-1.0||U_i[threadIdx.x]==-1.0)
////                U_i[threadIdx.x]=-1;
////            //add a partial sum upstream to our own
////            else
//            U_i[threadIdx.x] += U_i[threadIdx.x + offset];
////                U_i[threadIdx.x]=offset;
//        }
//        // wait until all threads in the block have updated their partial sums
//        __syncthreads();
//    }
//
//    __syncthreads();
//    // thread 0 writes the final result
//    if(threadIdx.x == 0)
//    {
//        U[blockIdx.x] = U_i[0];
//    }
//    __syncthreads();
    }
}


__global__ void CUDAExclVolRepulsionforce(double *coord, double *f, int *beadSet, double *krep, int *params){


    int nint = params[1];
    int n = params[0];
    double d, invDSquare, U;
    double a, b, c, e, F, AA, BB, CC, DD, EE, FF, GG, HH, JJ, invJJ;
    double ATG1, ATG2, ATG3, ATG4;
    double A1, A2, E1, E2, B1, B2, F1, F2, A11, A12, A13, A14;
    double E11, E12, E13, E14, B11, B12, B13, B14, F11, F12, F13, F14;
    double f1[3], f2[3], f3[3], f4[3];

    extern __shared__ double s[];
    double *c1 = s;
    double *c2 = &c1[3 * blockDim.x];
    double *c3 = &c2[3 * blockDim.x];
    double *c4 = &c3[3 * blockDim.x];
//    double *f1 = &c4[3 * blockDim.x];
//    double *f2 = &f1[3 * blockDim.x];
//    double *f3 = &f2[3 * blockDim.x];
//    double *f4 = &f3[3 * blockDim.x];

    const unsigned int thread_idx = (blockIdx.x * blockDim.x) + threadIdx.x;

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

//    }
//    __syncthreads();
//    if(thread_idx<nint) {

        //check if parallel
        if(areParallel(c1, c2, c3, c4, 3 * threadIdx.x))  {

            d = twoPointDistance(c1, c3, 3 * threadIdx.x);
            invDSquare =  1 / (d * d);
            U = krep[thread_idx] * invDSquare * invDSquare;

            double f0 = 4 * krep[thread_idx] * invDSquare * invDSquare * invDSquare;

            f1[0] = - f0 * (c3[0] - c1[0]);
            f1[1] = - f0 * (c3[1] - c1[1]);
            f1[2] = - f0 * (c3[2] - c1[2]);

            f2[0] = - f0 * (c4[0] - c2[0]);
            f2[1] = - f0 * (c4[1] - c2[1]);
            f2[2] = - f0 * (c4[2] - c2[2]);

            f3[0] = f0 * (c3[0] - c1[0]);
            f3[1] = f0 * (c3[1] - c1[1]);
            f3[2] = f0 * (c3[2] - c1[2]);

            f4[0] = f0 * (c4[0] - c2[0]);
            f4[1] = f0 * (c4[1] - c2[1]);
            f4[2] = f0 * (c4[2] - c2[2]);

            for(int i=0;i<3;i++) {
                atomicAdd(&f[3 * beadSet[n * thread_idx]  +i], f1[i]);
                atomicAdd(&f[3 * beadSet[n * thread_idx+1]  +i], f2[i]);
                atomicAdd(&f[3 * beadSet[n * thread_idx+2]  +i], f3[i]);
                atomicAdd(&f[3 * beadSet[n * thread_idx+3]  +i], f4[i]);
            }
        }

        else {

            //check if in same plane
            if(areInPlane(c1, c2, c3, c4, 3 * threadIdx.x)) {

                //slightly move point
                movePointOutOfPlane(c1, c2, c3, c4, 2, 0.01, 3 * threadIdx.x);
            }

            a = scalarProduct(c1, c2, c1, c2, 3 * threadIdx.x);
            b = scalarProduct(c3, c4, c3, c4, 3 * threadIdx.x);
            c = scalarProduct(c3, c1, c3, c1, 3 * threadIdx.x);
            d = scalarProduct(c1, c2, c3, c4, 3 * threadIdx.x);
            e = scalarProduct(c1, c2, c3, c1, 3 * threadIdx.x);
            F = scalarProduct(c3, c4, c3, c1, 3 * threadIdx.x);

            AA = sqrt(a*c - e*e);
            BB = sqrt(b*c - F*F);

            CC = d*e - a*F;
            DD = b*e - d*F;

            EE = sqrt( a*(b + c - 2*F) - (d - e)*(d - e) );
            FF = sqrt( b*(a + c + 2*e) - (d + F)*(d + F) );

            GG = d*d - a*b - CC;
            HH = CC + GG - DD;
            JJ = c*(GG + CC) + e*DD - F*CC;

            invJJ = 1/JJ;

            ATG1 = atan( (a + e)/AA) - atan(e/AA);
            ATG2 = atan((a + e - d)/EE) - atan((e - d)/EE);
            ATG3 = atan((F)/BB) - atan((F - b)/BB);
            ATG4 = atan((d + F)/FF) - atan((d + F - b)/FF);

            U = 0.5 * krep[thread_idx]/ JJ * ( CC/AA*ATG1 + GG/EE*ATG2 + DD/BB*ATG3 + HH/FF*ATG4);


            A1 = AA * AA / (AA * AA + e * e);
            A2 = AA * AA / (AA * AA + (a + e) * (a + e));

            E1 = EE * EE / (EE * EE + (a + e - d) * (a + e - d));
            E2 = EE * EE / (EE * EE + (e - d) * (e - d));

            B1 = BB * BB / (BB * BB + (F - b) * (F - b));;
            B2 = BB * BB / (BB * BB + F * F);

            F1 = FF * FF / (FF * FF + (d + F - b) * (d + F - b));
            F2 = FF * FF / (FF * FF + (d + F) * (d + F));

            A11 = ATG1 / AA;
            A12 = -((ATG1 * CC) / (AA * AA)) + (A1 * CC * e) / (AA * AA * AA) -
                  (A2 * CC * (a + e)) / (AA * AA * AA);
            A13 = -((A1 * CC) / (AA * AA)) + (A2 * CC) / (AA * AA);
            A14 = (A2 * CC) / (AA * AA);

            E11 = ATG2 / EE;
            E12 = (E2 * (-a + d - e) * GG) / (EE * EE * EE) + (E1 * (-d + e) * GG) / (EE * EE * EE) -
                  (ATG2 * GG) / (EE * EE);
            E13 = -((E1 * GG) / (EE * EE)) + (E2 * GG) / (EE * EE);
            E14 = (E2 * GG) / (EE * EE);

            B11 = ATG3 / BB;
            B12 = -((ATG3 * DD) / (BB * BB)) - (B2 * DD * F) / (BB * BB * BB) +
                  (B1 * DD * (-b + F)) / (BB * BB * BB);
            B13 = -((B1 * DD) / (BB * BB)) + (B2 * DD) / (BB * BB);
            B14 = (B1 * DD) / (BB * BB);

            F11 = ATG4 / FF;
            F12 = (F2 * (-d - F) * HH) / (FF * FF * FF) + (F1 * (-b + d + F) * HH) / (FF * FF * FF) -
                  (ATG4 * HH) / (FF * FF);
            F13 = -((F1 * HH) / (FF * FF)) + (F2 * HH) / (FF * FF);
            F14 = (F1 * HH) / (FF * FF);


            f1[0] =  - 0.5*invJJ*( (c2[3 * threadIdx.x] - c1[3 * threadIdx.x] ) *( A13 + E13 + B11*b - F11*b + A11*d - E11*d - 2*U*b*e - (A12*e)/AA + (E12*(d - e))/EE + 2*U*d*F - 2*U*(b*e - d*F) + (F12*b)/FF - 2*(A14 + E14 - E11*b - F11*b + 2*U*b*c + (A12*c)/(2*AA) + (E12*(b + c - 2*F))/(2*EE) - A11*F + E11*F - 2*U*F*F + (F12*b)/(2*FF)) ) + (c4[3 * threadIdx.x] - c3[3 * threadIdx.x] ) *(B13 + E13 - A11*a + E11*a - B11*d - 2*E11*d - F11*d + 4*U*c*d - A11*e + E11*e + 2*U*d*e - (E12*a)/EE + (E12*(d - e))/EE + B11*F - F11*F - 2*U*a*F - 4*U*e*F + 2*U*(d*e - a*F) - (B12*F)/BB) +  (c1[3 * threadIdx.x] - c3[3 * threadIdx.x] )* (-A13 - E13 - B11*b + F11*b - A11*d + E11*d + 2*U*b*e + (A12*e)/AA - (E12*(d - e))/EE - 2*U*d*F + 2*U*(b*e - d*F) - (F12*b)/FF + 2*(-2*U*((-a)*b + d*d) + (A12*a)/(2*AA) + (E12*a)/(2*EE) +(B12*b)/(2*BB) + (F12*b)/(2*FF))) );

            f1[1] =  - 0.5*invJJ*( (c2[3 * threadIdx.x + 1] - c1[3 * threadIdx.x + 1] ) *( A13 + E13 + B11*b - F11*b + A11*d - E11*d - 2*U*b*e - (A12*e)/AA + (E12*(d - e))/EE + 2*U*d*F - 2*U*(b*e - d*F) + (F12*b)/FF - 2*(A14 + E14 - E11*b - F11*b + 2*U*b*c + (A12*c)/(2*AA) + (E12*(b + c - 2*F))/(2*EE) - A11*F + E11*F - 2*U*F*F + (F12*b)/(2*FF)) ) + (c4[3 * threadIdx.x + 1] - c3[3 * threadIdx.x + 1] ) *(B13 + E13 - A11*a + E11*a - B11*d - 2*E11*d - F11*d + 4*U*c*d - A11*e + E11*e + 2*U*d*e - (E12*a)/EE + (E12*(d - e))/EE + B11*F - F11*F - 2*U*a*F - 4*U*e*F + 2*U*(d*e - a*F) - (B12*F)/BB) +  (c1[3 * threadIdx.x + 1] - c3[3 * threadIdx.x + 1] )* (-A13 - E13 - B11*b + F11*b - A11*d + E11*d + 2*U*b*e + (A12*e)/AA - (E12*(d - e))/EE - 2*U*d*F + 2*U*(b*e - d*F) - (F12*b)/FF + 2*(-2*U*((-a)*b + d*d) + (A12*a)/(2*AA) + (E12*a)/(2*EE) +(B12*b)/(2*BB) + (F12*b)/(2*FF))) );

            f1[2] =  - 0.5*invJJ*( (c2[3 * threadIdx.x + 2] - c1[3 * threadIdx.x + 2] ) *( A13 + E13 + B11*b - F11*b + A11*d - E11*d - 2*U*b*e - (A12*e)/AA + (E12*(d - e))/EE + 2*U*d*F - 2*U*(b*e - d*F) + (F12*b)/FF - 2*(A14 + E14 - E11*b - F11*b + 2*U*b*c + (A12*c)/(2*AA) + (E12*(b + c - 2*F))/(2*EE) - A11*F + E11*F - 2*U*F*F + (F12*b)/(2*FF)) ) + (c4[3 * threadIdx.x + 2] - c3[3 * threadIdx.x + 2] ) *(B13 + E13 - A11*a + E11*a - B11*d - 2*E11*d - F11*d + 4*U*c*d - A11*e + E11*e + 2*U*d*e - (E12*a)/EE + (E12*(d - e))/EE + B11*F - F11*F - 2*U*a*F - 4*U*e*F + 2*U*(d*e - a*F) - (B12*F)/BB) +  (c1[3 * threadIdx.x + 2] - c3[3 * threadIdx.x + 2] )* (-A13 - E13 - B11*b + F11*b - A11*d + E11*d + 2*U*b*e + (A12*e)/AA - (E12*(d - e))/EE - 2*U*d*F + 2*U*(b*e - d*F) - (F12*b)/FF + 2*(-2*U*((-a)*b + d*d) + (A12*a)/(2*AA) + (E12*a)/(2*EE) +(B12*b)/(2*BB) + (F12*b)/(2*FF))) );


            f2[0] =  - invJJ*( (c2[3 * threadIdx.x] - c1[3 * threadIdx.x] )*( A14+E14-E11*b-F11*b+2*U*b*c+(A12*c)/(2*AA)+(E12*(b+c-2*F))/(2*EE)-A11*F+E11*F-2*U*F*F+(F12*b)/(2*FF) ) + 0.5*(c4[3 * threadIdx.x] - c3[3 * threadIdx.x])*(-E13 + F13 + 2*E11*d + 2*F11*d - 4*U*c*d + A11*e - E11*e - (E12*(d - e))/EE - B11*F + F11*F + 4*U*e*F - (F12*(d + F))/FF)  + 0.5*(c1[3 * threadIdx.x] - c3[3 * threadIdx.x] )* (A13 + E13 + B11*b - F11*b + A11*d - E11*d - 2*U*b*e - (A12*e)/AA + (E12*(d - e))/EE + 2*U*d*F - 2*U*(b*e - d*F) + (F12*b)/FF) );

            f2[1] = - invJJ*( (c2[3 * threadIdx.x + 1] - c1[3 * threadIdx.x + 1] )*( A14+E14-E11*b-F11*b+2*U*b*c+(A12*c)/(2*AA)+(E12*(b+c-2*F))/(2*EE)-A11*F+E11*F-2*U*F*F+(F12*b)/(2*FF) ) + 0.5*(c4[3 * threadIdx.x + 1] - c3[3 * threadIdx.x + 1])*(-E13 + F13 + 2*E11*d + 2*F11*d - 4*U*c*d + A11*e - E11*e - (E12*(d - e))/EE - B11*F + F11*F + 4*U*e*F - (F12*(d + F))/FF)  + 0.5*(c1[3 * threadIdx.x + 1] - c3[3 * threadIdx.x + 1] )* (A13 + E13 + B11*b - F11*b + A11*d - E11*d - 2*U*b*e - (A12*e)/AA + (E12*(d - e))/EE + 2*U*d*F - 2*U*(b*e - d*F) + (F12*b)/FF) );

            f2[2] = - invJJ*( (c2[3 * threadIdx.x + 2] - c1[3 * threadIdx.x + 2] )*( A14+E14-E11*b-F11*b+2*U*b*c+(A12*c)/(2*AA)+(E12*(b+c-2*F))/(2*EE)-A11*F+E11*F-2*U*F*F+(F12*b)/(2*FF) ) + 0.5*(c4[3 * threadIdx.x + 2] - c3[3 * threadIdx.x + 2])*(-E13 + F13 + 2*E11*d + 2*F11*d - 4*U*c*d + A11*e - E11*e - (E12*(d - e))/EE - B11*F + F11*F + 4*U*e*F - (F12*(d + F))/FF)  + 0.5*(c1[3 * threadIdx.x + 2] - c3[3 * threadIdx.x + 2] )* (A13 + E13 + B11*b - F11*b + A11*d - E11*d - 2*U*b*e - (A12*e)/AA + (E12*(d - e))/EE + 2*U*d*F - 2*U*(b*e - d*F) + (F12*b)/FF) );

            f3[0] =  - 0.5*invJJ*( (c2[3 * threadIdx.x] - c1[3 * threadIdx.x] )*(-A13 - F13 - B11*b + F11*b - A11*d - E11*d - 2*F11*d + 4*U*c*d - A11*e + E11*e + 2*U*b*e + (A12*e)/AA + B11*F - F11*F - 2*U*d*F - 4*U*e*F + 2*U*(b*e - d*F) - (F12*b)/FF + (F12*(d + F))/FF) + (c4[3 * threadIdx.x] - c3[3 * threadIdx.x] )*(-B13 - F13 + A11*a - E11*a + B11*d - F11*d - 2*U*d*e + (E12*a)/EE + 2*U*a*F - 2*U*(d*e - a*F) + (B12*F)/BB + (F12*(d + F))/FF - 2*(B14 + F14 - E11*a - F11*a + 2*U*a*c + B11*e - F11*e - 2*U*e*e + (E12*a)/(2*EE) + (B12*c)/(2*BB) + (F12*(a + c + 2*e))/(2*FF))) + (c1[3 * threadIdx.x] - c3[3 * threadIdx.x] ) * (-B13 - F13 + A11*a - E11*a + B11*d - F11*d - 2*U*d*e + (E12*a)/EE + 2*U*a*F - 2*U*(d*e - a*F) + (B12*F)/BB + (F12*(d + F))/FF - 2*(-2*U*((-a)*b + d*d) + (A12*a)/(2*AA) + (E12* a)/(2*EE) + (B12*b)/(2*BB) + (F12*b)/(2*FF))) );

            f3[1] =  - 0.5*invJJ*( (c2[3 * threadIdx.x + 1] - c1[3 * threadIdx.x + 1] )*(-A13 - F13 - B11*b + F11*b -
                    A11*d - E11*d - 2*F11*d + 4*U*c*d - A11*e + E11*e + 2*U*b*e + (A12*e)/AA + B11*F - F11*F - 2*U*d*F - 4*U*e*F + 2*U*(b*e - d*F) - (F12*b)/FF + (F12*(d + F))/FF) + (c4[3 * threadIdx.x + 1] - c3[3 * threadIdx.x + 1] )*(-B13 - F13 + A11*a - E11*a + B11*d - F11*d - 2*U*d*e + (E12*a)/EE + 2*U*a*F - 2*U*(d*e - a*F) + (B12*F)/BB + (F12*(d + F))/FF - 2*(B14 + F14 - E11*a - F11*a + 2*U*a*c + B11*e - F11*e - 2*U*e*e + (E12*a)/(2*EE) + (B12*c)/(2*BB) + (F12*(a + c + 2*e))/(2*FF))) + (c1[3 * threadIdx.x + 1] - c3[3 * threadIdx.x + 1] ) * (-B13 - F13 + A11*a - E11*a + B11*d - F11*d - 2*U*d*e + (E12*a)/EE + 2*U*a*F - 2*U*(d*e - a*F) + (B12*F)/BB + (F12*(d + F))/FF - 2*(-2*U*((-a)*b + d*d) + (A12*a)/(2*AA) + (E12* a)/(2*EE) + (B12*b)/(2*BB) + (F12*b)/(2*FF))) ) ;

            f3[2] =  - 0.5*invJJ*( (c2[3 * threadIdx.x + 2] - c1[3 * threadIdx.x + 2] )*(-A13 - F13 - B11*b + F11*b - A11*d - E11*d - 2*F11*d + 4*U*c*d - A11*e + E11*e + 2*U*b*e + (A12*e)/AA + B11*F - F11*F - 2*U*d*F - 4*U*e*F + 2*U*(b*e - d*F) - (F12*b)/FF + (F12*(d + F))/FF) + (c4[3 * threadIdx.x + 2] - c3[3 * threadIdx.x + 2] )*(-B13 - F13 + A11*a - E11*a + B11*d - F11*d - 2*U*d*e + (E12*a)/EE + 2*U*a*F - 2*U*(d*e - a*F) + (B12*F)/BB + (F12*(d + F))/FF - 2*(B14 + F14 - E11*a - F11*a + 2*U*a*c + B11*e - F11*e - 2*U*e*e + (E12*a)/(2*EE) + (B12*c)/(2*BB) + (F12*(a + c + 2*e))/(2*FF))) + (c1[3 * threadIdx.x + 2] - c3[3 * threadIdx.x + 2] ) * (-B13 - F13 + A11*a - E11*a + B11*d - F11*d - 2*U*d*e + (E12*a)/EE + 2*U*a*F - 2*U*(d*e - a*F) + (B12*F)/BB + (F12*(d + F))/FF - 2*(-2*U*((-a)*b + d*d) + (A12*a)/(2*AA) + (E12* a)/(2*EE) + (B12*b)/(2*BB) + (F12*b)/(2*FF))) );


            f4[0] =  - invJJ*( 0.5*(c2[3 * threadIdx.x] - c1[3 * threadIdx.x] )*( -E13 + F13 + 2*E11*d + 2*F11*d - 4*U*c*d + A11*e - E11*e - (E12*(d - e))/EE - B11*F + F11*F +4*U*e*F - (F12*(d + F))/FF ) + (c4[3 * threadIdx.x] - c3[3 * threadIdx.x])*(B14 + F14 - E11*a - F11*a + 2*U*a*c + B11*e - F11*e - 2*U*e*e + (E12*a)/(2*EE) + (B12*c)/(2*BB) + (F12*(a + c + 2*e))/(2*FF))  + 0.5*(c1[3 * threadIdx.x] - c3[3 * threadIdx.x] )* (B13 + F13 - A11*a + E11*a - B11*d + F11*d + 2*U*d*e - (E12*a)/EE - 2*U*a*F + 2*U*(d*e - a*F) - (B12*F)/BB - (F12*(d + F))/FF) )  ;

            f4[1] =  - invJJ*( 0.5*(c2[3 * threadIdx.x + 1] - c1[3 * threadIdx.x + 1] )*( -E13 + F13 + 2*E11*d + 2*F11*d - 4*U*c*d + A11*e - E11*e - (E12*(d - e))/EE - B11*F + F11*F +4*U*e*F - (F12*(d + F))/FF ) + (c4[3 * threadIdx.x + 1] - c3[3 * threadIdx.x + 1])*(B14 + F14 - E11*a - F11*a + 2*U*a*c + B11*e - F11*e - 2*U*e*e + (E12*a)/(2*EE) + (B12*c)/(2*BB) + (F12*(a + c + 2*e))/(2*FF))  + 0.5*(c1[3 * threadIdx.x + 1] - c3[3 * threadIdx.x + 1] )* (B13 + F13 - A11*a + E11*a - B11*d + F11*d + 2*U*d*e - (E12*a)/EE - 2*U*a*F + 2*U*(d*e - a*F) - (B12*F)/BB - (F12*(d + F))/FF) ) ;

            f4[2] =  - invJJ*( 0.5*(c2[3 * threadIdx.x + 2] - c1[3 * threadIdx.x + 2] )*( -E13 + F13 + 2*E11*d + 2*F11*d - 4*U*c*d + A11*e - E11*e - (E12*(d - e))/EE - B11*F + F11*F +4*U*e*F - (F12*(d + F))/FF ) + (c4[3 * threadIdx.x + 2] - c3[3 * threadIdx.x + 2])*(B14 + F14 - E11*a - F11*a + 2*U*a*c + B11*e - F11*e - 2*U*e*e + (E12*a)/(2*EE) + (B12*c)/(2*BB) + (F12*(a + c + 2*e))/(2*FF))  + 0.5*(c1[3 * threadIdx.x + 2] - c3[3 * threadIdx.x + 2] )* (B13 + F13 - A11*a + E11*a - B11*d + F11*d + 2*U*d*e - (E12*a)/EE - 2*U*a*F + 2*U*(d*e - a*F) - (B12*F)/BB - (F12*(d + F))/FF) ) ;


            for(int i=0;i<3;i++) {
                atomicAdd(&f[3 * beadSet[n * thread_idx]  +i], f1[i]);
                atomicAdd(&f[3 * beadSet[n * thread_idx+1]  +i], f2[i]);
                atomicAdd(&f[3 * beadSet[n * thread_idx+2]  +i], f3[i]);
                atomicAdd(&f[3 * beadSet[n * thread_idx+3]  +i], f4[i]);
            }

        }
//        for(int i=0;i<3;i++) {
//            f1c[3 * thread_idx + i] =  f1[i];
//            f2c[3 * thread_idx + i] = f2[i];
//            f3c[3 * thread_idx + i] =  f3[i];
//            f4c[3 * thread_idx + i] =  f4[i];
//        }

//            f1c[3 * thread_idx]  =  f1[0];
//        f1c[3 * thread_idx + 1]  = f1[1];
//        f1c[3 * thread_idx + 2]  =  f1[2];
//        f2c[3 * thread_idx]  =  d;
//        f2c[3 * thread_idx + 1]  = e;
//        f2c[3 * thread_idx + 2]  =  F;
//        f3c[3 * thread_idx]  =  AA;
//        f3c[3 * thread_idx + 1]  =  BB;
//        f3c[3 * thread_idx + 2]  =  A1;
//        f4c[3 * thread_idx]  =  A2;
//        f4c[3 * thread_idx + 1]  = 0.0;
//        f4c[3 * thread_idx + 2]  =  0.0;

//        f5c[44 * thread_idx]  =  a;
//        f5c[44 * thread_idx+1]  =  b;
//        f5c[44 * thread_idx+2]  =  c;
//        f5c[44 * thread_idx+3]  =  d;
//        f5c[44 * thread_idx+4]  =  e;
//        f5c[44 * thread_idx+5]  =  F;
//        f5c[44 * thread_idx+6]  =  AA;
//        f5c[44 * thread_idx+7]  =  BB;
//        f5c[44 * thread_idx+8]  =  CC;
//        f5c[44 * thread_idx+9]  =  DD;
//        f5c[44 * thread_idx+10]  =  EE;
//        f5c[44 * thread_idx+11]  =  FF;
//        f5c[44 * thread_idx+12]  =  GG;
//        f5c[44 * thread_idx+13]  =  HH;
//        f5c[44 * thread_idx+14]  =  JJ;
//        f5c[44 * thread_idx+15]  =  ATG1;
//        f5c[44 * thread_idx+16]  =  ATG2;
//        f5c[44 * thread_idx+17]  =  ATG3;
//        f5c[44 * thread_idx+18]  =  ATG4;
//        f5c[44 * thread_idx+19]  =  U;
//        f5c[44 * thread_idx+20]  =  A1;
//        f5c[44 * thread_idx+21]  =  A2;
//        f5c[44 * thread_idx+22]  =  E1;
//        f5c[44 * thread_idx+23]  =  E2;
//        f5c[44 * thread_idx+24]  =  B1;
//        f5c[44 * thread_idx+25]  =  B2;
//        f5c[44 * thread_idx+26]  =  F1;
//        f5c[44 * thread_idx+27]  =  F2;
//        f5c[44 * thread_idx+28]  =  A11;
//        f5c[44 * thread_idx+29]  =  A12;
//        f5c[44 * thread_idx+30]  =  A13;
//        f5c[44 * thread_idx+31]  =  A14;
//        f5c[44 * thread_idx+32]  =  E11;
//        f5c[44 * thread_idx+33]  =  E12;
//        f5c[44 * thread_idx+34]  =  E13;
//        f5c[44 * thread_idx+35]  =  E14;
//        f5c[44 * thread_idx+36]  =  B11;
//        f5c[44 * thread_idx+37]  =  B12;
//        f5c[44 * thread_idx+38]  =  B13;
//        f5c[44 * thread_idx+39]  =  B14;
//        f5c[44 * thread_idx+40]  =  F11;
//        f5c[44 * thread_idx+41]  =  F12;
//        f5c[44 * thread_idx+42]  =  F13;
//        f5c[44 * thread_idx+43]  =  F14;


    }

}



#endif
#endif //CUDA_VEC_CYLINDEREXCLVOLREPULSIONCUDA_H
