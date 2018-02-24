//
// Created by aravind on 11/21/17.
//

#ifndef CUDA_VEC_BOUNDARYCYLINDERREPULSIONEXPCUDA_H
#define CUDA_VEC_BOUNDARYCYLINDERREPULSIONEXPCUDA_H

#include "MathFunctions.h"
#include "SysParams.h"
#include <limits>
#include <assert.h>
#include <cuda.h>
#include <cuda_runtime.h>
using namespace mathfunc;
__global__ void BoundaryCylinderRepulsionExpenergy(double* coord, double* f, int* beadSet, double* krep, double* slen,
                                                   int* nintvec, double* beListplane, int* params, double* U_i, double *z,
                                                   int* culpritID, char* culpritFF, char* culpritinteraction, char*
                                                   FF, char* interaction){
    if(z[0] == 0.0) {
        extern __shared__ double s[];
        double *c1 = s;
        int nint = params[1];
        double R, r;
        int n = params[0];
        const unsigned int thread_idx = (blockIdx.x * blockDim.x) + threadIdx.x;
        double plane[4];
        if (thread_idx < nint) {
            U_i[thread_idx] = 0.0;
            for (auto i = 0; i < 3; i++) {
                U_i[thread_idx] = 0.0;
                c1[3 * threadIdx.x + i] = coord[3 * beadSet[n * thread_idx] + i];
            }

        }
        __syncthreads();
        if (thread_idx < nint) {
            //get the plane equation.
            if (thread_idx < nintvec[0]) {
                for (auto i = 0; i < 3; i++)
                    plane[i] = beListplane[i];
            } else if (thread_idx >= nintvec[0] || thread_idx < nintvec[1]) {
                for (auto i = 0; i < 4; i++)
                    plane[i] = beListplane[4 + i];
            } else if (thread_idx >= nintvec[1] || thread_idx < nintvec[2]) {
                for (auto i = 0; i < 4; i++)
                    plane[i] = beListplane[8 + i];
            } else if (thread_idx >= nintvec[2] || thread_idx < nintvec[3]) {
                for (auto i = 0; i < 4; i++)
                    plane[i] = beListplane[12 + i];
            } else if (thread_idx >= nintvec[3] || thread_idx < nintvec[4]) {
                for (auto i = 0; i < 4; i++)
                    plane[i] = beListplane[16 + i];
            } else if (thread_idx >= nintvec[4] || thread_idx < nintvec[5]) {
                for (auto i = 0; i < 4; i++)
                    plane[i] = beListplane[20 + i];
            }
            //get distance from plane
            r = getdistancefromplane(c1, plane, 3 * threadIdx.x);
            R = -r / slen[thread_idx];
            U_i[thread_idx] = krep[thread_idx] * exp(R);
//        printf("%d %f %f\n", thread_idx,krep[thread_idx], R);
            if (fabs(U_i[thread_idx]) == __longlong_as_double(0x7ff0000000000000) //infinity
                || U_i[thread_idx] != U_i[thread_idx] || U_i[thread_idx] < -1.0) {
                //set culprit and exit
                U_i[thread_idx] = -1.0;
                culpritID[0] = beadSet[n * thread_idx];//set Cylinder info.
                //set boundary element information
                if (thread_idx < nintvec[0]) culpritID[1] = 0;
                else if (thread_idx >= nintvec[0] || thread_idx < nintvec[1]) culpritID[1] = 1;
                else if (thread_idx >= nintvec[1] || thread_idx < nintvec[2]) culpritID[1] = 2;
                else if (thread_idx >= nintvec[2] || thread_idx < nintvec[3]) culpritID[1] = 3;
                else if (thread_idx >= nintvec[3] || thread_idx < nintvec[4]) culpritID[1] = 4;
                else if (thread_idx >= nintvec[4] || thread_idx < nintvec[5]) culpritID[1] = 5;
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
                printf("%d %f %f %f\n", thread_idx, r, R, U_i[thread_idx]);
                assert(0);
                __syncthreads();
            }
        }
    }
}

__global__ void BoundaryCylinderRepulsionExpenergyz(double* coord, double* f, int* beadSet, double* krep, double* slen,
                                                   int* nintvec, double* beListplane, int* params, double* U_i,
                                                   double *z, int* culpritID, char* culpritFF, char* culpritinteraction,
                                                   char* FF, char* interaction){
    if(z[0] != 0.0) {
        extern __shared__ double s[];
        double *c1 = s;
        double *f1 = &c1[3 * blockDim.x];
        int nint = params[1];
        double R, r;
        int n = params[0];
        const unsigned int thread_idx = (blockIdx.x * blockDim.x) + threadIdx.x;
        double plane[4];
        if (thread_idx < nint) {
            U_i[thread_idx] = 0.0;
            for (auto i = 0; i < 3; i++) {
                U_i[thread_idx] = 0.0;
                c1[3 * threadIdx.x + i] = coord[3 * beadSet[n * thread_idx] + i];
                f1[3 * threadIdx.x + i] = f[3 * beadSet[n * thread_idx] + i];
            }
        }
        __syncthreads();
        if (thread_idx < nint) {
            //get the plane equation.
            if (thread_idx < nintvec[0]) {
                for (auto i = 0; i < 3; i++)
                    plane[i] = beListplane[i];
            } else if (thread_idx >= nintvec[0] || thread_idx < nintvec[1]) {
                for (auto i = 0; i < 4; i++)
                    plane[i] = beListplane[4 + i];
            } else if (thread_idx >= nintvec[1] || thread_idx < nintvec[2]) {
                for (auto i = 0; i < 4; i++)
                    plane[i] = beListplane[8 + i];
            } else if (thread_idx >= nintvec[2] || thread_idx < nintvec[3]) {
                for (auto i = 0; i < 4; i++)
                    plane[i] = beListplane[12 + i];
            } else if (thread_idx >= nintvec[3] || thread_idx < nintvec[4]) {
                for (auto i = 0; i < 4; i++)
                    plane[i] = beListplane[16 + i];
            } else if (thread_idx >= nintvec[4] || thread_idx < nintvec[5]) {
                for (auto i = 0; i < 4; i++)
                    plane[i] = beListplane[20 + i];
            }
            //get distance from plane
            r = getstretcheddistancefromplane(c1, f1, plane, z[0], 3 * threadIdx.x);
            R = -r / slen[thread_idx];
            U_i[thread_idx] = krep[thread_idx] * exp(R);
//        printf("Z %d %f %f\n", thread_idx,krep[thread_idx], R);
//        printf("%d %f %f %f\n", thread_idx, r, R, U_i[thread_idx]);
            if (fabs(U_i[thread_idx]) == __longlong_as_double(0x7ff0000000000000) //infinity
                || U_i[thread_idx] != U_i[thread_idx] || U_i[thread_idx] < -1.0) {
                //set culprit and exit
                U_i[thread_idx] = -1.0;
                culpritID[0] = thread_idx;
                if (thread_idx < nintvec[0]) culpritID[1] = 0;
                else if (thread_idx >= nintvec[0] || thread_idx < nintvec[1]) culpritID[1] = 1;
                else if (thread_idx >= nintvec[1] || thread_idx < nintvec[2]) culpritID[1] = 2;
                else if (thread_idx >= nintvec[2] || thread_idx < nintvec[3]) culpritID[1] = 3;
                else if (thread_idx >= nintvec[3] || thread_idx < nintvec[4]) culpritID[1] = 4;
                else if (thread_idx >= nintvec[4] || thread_idx < nintvec[5]) culpritID[1] = 5;
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
                printf("%d %f %f %f\n", thread_idx, r, R, U_i[thread_idx]);
                assert(0);
                __syncthreads();
            }
        }
    }
}

__global__ void BoundaryCylinderRepulsionExpforces(double* coord, double* f, int* beadSet, double* krep, double* slen,
                                                   int* nintvec, double* beListplane, int* params){
    extern __shared__ double s[];
    double *c1 = s;
    int nint = params[1];
    double R, r, norm[3], f0;
    int n = params[0];
    const unsigned int thread_idx = (blockIdx.x * blockDim.x) + threadIdx.x;
    double plane[4];
    if(thread_idx<nint) {
        for(auto i=0;i<3;i++){
            c1[3 * threadIdx.x + i] = coord[3 * beadSet[n * thread_idx] + i];
        }

    }
    __syncthreads();
    if(thread_idx<nint) {
        //get the plane equation.
        if (thread_idx < nintvec[0]) {
            for (auto i = 0; i<3; i++)
                plane[i] = beListplane [i];
        }
        else if(thread_idx >= nintvec[0] || thread_idx < nintvec[1]){
            for (auto i = 0; i<4; i++)
                plane[i] = beListplane [4 + i];
        }
        else if(thread_idx >= nintvec[1] || thread_idx < nintvec[2]){
            for (auto i = 0; i<4; i++)
                plane[i] = beListplane [8 + i];
        }
        else if(thread_idx >= nintvec[2] || thread_idx < nintvec[3]){
            for (auto i = 0; i<4; i++)
                plane[i] = beListplane [12 + i];
        }
        else if(thread_idx >= nintvec[3] || thread_idx < nintvec[4]){
            for (auto i = 0; i<4; i++)
                plane[i] = beListplane [16 + i];
        }
        else if(thread_idx >= nintvec[4] || thread_idx < nintvec[5]){
            for (auto i = 0; i<4; i++)
                plane[i] = beListplane [20 + i];
        }
        //get distance from plane
        r = getdistancefromplane(c1, plane, 3 * threadIdx.x);
        for(auto i = 0; i < 3; i++)
            norm[i] = plane[i];
        R = -r / slen[thread_idx];
        f0 = krep[thread_idx] * exp(R);
        for (int i = 0; i < 3; i++) {
            atomicAdd(&f[3 * beadSet[n * thread_idx] + i], f0*norm[i]);
        }
    }
}
//__global__ void BoundaryCylinderRepulsionadd(double *force, double *forcecopy, int *nint) {
//    const unsigned int thread_idx = (blockIdx.x * blockDim.x) + threadIdx.x;
//
//    if (thread_idx < nint[0]) {
//        for (auto i = 0; i < 3; i++) {
//            force[3 * thread_idx + i] =  force[3 * thread_idx + i] + forcecopy[3 * thread_idx + i];
//        }
//    }
//}

#endif //CUDA_VEC_BOUNDARYCYLINDERREPULSIONCUDA_H
