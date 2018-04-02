
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

#include <cmath>
#include <vector>
#include <math.h>

#include "MathFunctions.h"
#include "Rand.h"

namespace mathfunc {

    tuple<vector<double>, vector<double>> branchProjection(const vector<double>& n,
                                                           const vector<double>& p,
                                                           double l, double m, double theta){
        //get random permutation from p
        vector<double> r = {p[0] + Rand::randDouble(-1, 1),
                            p[1] + Rand::randDouble(-1, 1),
                            p[2] + Rand::randDouble(-1, 1)};

        //construct vector z which is r-p
        auto z = twoPointDirection(p, r);

        //construct u and v, which creates an orthogonal set n, u, v
        auto u = crossProduct(n, z);
        auto v = crossProduct(n, u);

        normalize(u); normalize(v);

        //find random point on circle defining the branching point
        double thetaRandom = Rand::randDouble(0, 2*M_PI);
        vector<double> bp1;
        bp1.push_back(p[0] + l * (u[0] * cos(thetaRandom) + v[0] * sin(thetaRandom)));
        bp1.push_back(p[1] + l * (u[1] * cos(thetaRandom) + v[1] * sin(thetaRandom)));
        bp1.push_back(p[2] + l * (u[2] * cos(thetaRandom) + v[2] * sin(thetaRandom)));

        //now find the second point
        vector<double> newP;
        double dist = m * cos(theta);
        newP.push_back(p[0] + n[0] * dist);
        newP.push_back(p[1] + n[1] * dist);
        newP.push_back(p[2] + n[2] * dist);
        double newL = (l + m * sin(theta));

        vector<double> bp2;
        bp2.push_back(newP[0] + newL * (u[0] * cos(thetaRandom) + v[0] * sin(thetaRandom)));
        bp2.push_back(newP[1] + newL * (u[1] * cos(thetaRandom) + v[1] * sin(thetaRandom)));
        bp2.push_back(newP[2] + newL * (u[2] * cos(thetaRandom) + v[2] * sin(thetaRandom)));

        //get direction
        auto direction = twoPointDirection(bp1, bp2);

        return tuple<vector<double>, vector<double>>(direction, bp1);
    }
     __global__ void addvector(double *U, int *params, double *U_sum, double *U_tot){
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
        printf("ADD1 %f\n", sum);

    }

//    __global__ void addvectorred(double *U, int *params, double *U_sum, double *U_tot){
//        extern __shared__ double s[];
//        double *c1 = s;
//        int start = 0;
//        int end = params[1];
//        int factor = params[1]/blockDim.x;
//        if(factor > 0) {
//            if (threadIdx.x > 0)
//                start = threadIdx.x * factor;
//            if (threadIdx.x < blockDim.x - 1)
//                end = (threadIdx.x + 1) * factor;
//            c1[threadIdx.x] = 0.0;
//            for (auto i = start; i < end; i++) {
//                if (c1[threadIdx.x] != -1.0 || U[i] != -1.0)
//                    c1[threadIdx.x] += U[i];
//                else
//                    c1[threadIdx.x] = -1.0;
//            }
////    printf("%d \n", params[0]);
//            __syncthreads();
//
//            if (threadIdx.x == 0) {
//                U_sum[0] = 0.0;
//                for (auto i = 0; i < blockDim.x; i++) {
//                    if (c1[i] != -1.0 && U_sum[0] != -1.0)
//                        U_sum[0] += c1[i];
//                    else
//                        U_sum[0] = -1.0;
//                }
//            if(U_sum[0] == -1.0){
//                auto val = -U_tot[0]-1.0;
//                atomicAdd(&U_tot[0], val);
//            }
//            else
//                atomicAdd(&U_tot[0], U_sum[0]);
////                printf("ADD2 %f %f \n", U_tot[0], U_sum[0]);
//            }
//        }
//        else if(threadIdx.x == 0) {
//            U_sum[0] = 0.0;
//            double sum = 0.0;
//            for (auto i = 0; i < params[1]; i++) {
//                if (U[i] == -1.0 && sum != -1.0) {
//                    U_sum[0] = -1.0;
//                    U_tot[0] = -1.0;
//                    sum = -1.0;
//                    break;
//                } else
//                    sum += U[i];
//            }
//            U_sum[0] = sum;
//            atomicAdd(&U_tot[0], sum);
////            printf("ADD3 %f %f \n", U_tot[0], U_sum[0]);
//        }
//    }

    __global__ void addvectorred2(double *g_idata, int *num, double *U_sum, double *g_odata)
    {
        extern __shared__ double sdata[];
        unsigned int tid = threadIdx.x;
        auto blockSize = blockDim.x;
        unsigned int i = blockIdx.x*(blockSize*2) + tid;
        unsigned int gridSize = blockSize*2*gridDim.x;
        int n = num[1];
        sdata[tid] = 0;
        while (i < n) {
//            printf("%d %d %d\n", tid, i, blockSize);
            if(g_idata[i] == -1.0 || g_idata[i+blockSize] == -1.0)
            {printf("CUDA addition of energies. Energy is infinite\n");assert(0);}
            sdata[tid] += g_idata[i] + g_idata[i+blockSize];
            i += gridSize;
        }
        __syncthreads();
        if(blockSize >=2048 ) {printf("Cannot handle blocks with threads larger than 2048\n");assert(0);}
        if (blockSize >= 1024) { if (tid < 512) { sdata[tid] += sdata[tid + 512]; } __syncthreads(); }
        if (blockSize >= 512) { if (tid < 256) { sdata[tid] += sdata[tid + 256]; } __syncthreads(); }
        if (blockSize >= 256) { if (tid < 128) { sdata[tid] += sdata[tid + 128]; } __syncthreads(); }
        if (blockSize >= 128) { if (tid < 64) { sdata[tid] += sdata[tid + 64]; } __syncthreads(); }
        if (tid < 32) {
            if (blockSize >= 64) sdata[tid] += sdata[tid + 32];
            if (blockSize >= 32) sdata[tid] += sdata[tid + 16];
            if (blockSize >= 16) sdata[tid] += sdata[tid + 8];
            if (blockSize >= 8) sdata[tid] += sdata[tid + 4];
            if (blockSize >= 4) sdata[tid] += sdata[tid + 2];
            if (blockSize >= 2) sdata[tid] += sdata[tid + 1];
        }
        if (tid == 0) {
            atomicAdd(&g_odata[0], sdata[0]);
            atomicAdd(&U_sum[0], sdata[0]);
//        printf("addv2 %f %f %f\n", sdata[0], sdata[1], sdata[2]);
        }
    }
    __global__ void resetintvariableCUDA(int *variable){
        variable[0] = 0;
    }

    __global__ void resetdoublevariableCUDA(double *variable){
        variable[0] = 0.0;
    }
//    __global__ void addvector(double *U, int *params, double *U_sum, double *U_tot, int *culpritID, char* culpritFF,
//                              char* culpritinteraction, char* FF, char* interaction){
//        U_sum[0] = 0.0;
//        double sum = 0.0;
//        for(auto i=0;i<params[1];i++){
//            if(U[i] == -1.0 && sum != -1.0){
//                U_sum[0] = -1.0;
//                U_tot[0] = -1.0;
//                sum = -1.0;
//            }
//            else if(sum != -1.0)
//                sum  += U[i];
//        }
//        U_sum[0] = sum;
////        printf("sum %f\n",sum);
//        if(sum != -1.0)
//            atomicAdd(&U_tot[0], sum);
//        else{
//            assert(0);
//        }
//    }
}