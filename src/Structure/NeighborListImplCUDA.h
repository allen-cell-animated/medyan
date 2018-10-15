//
// Created by aravind on 2/27/18.
//

#ifndef CUDA_VEC_NEIGHBORLISTIMPLCUDA_H
#define CUDA_VEC_NEIGHBORLISTIMPLCUDA_H
#ifdef CUDAACCL_TEST
#include <cuda.h>
#include <cuda_runtime.h>
#include <cuda_runtime_api.h>
#include <cmath>
#include "math.h"
#include "utility.h"
#include "assert.h"
#include "MathFunctions.h"
using namespace mathfunc;
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
//__global__ void testfunction(int *a){
//    int b =10.0;
//    b = atomicAdd(&a[0], 1);
//    printf("%d %d\n",b,a[0]);
//}

__global__ void CylinderCylinderNLCUDA(double *coord_com, int *beadSet, int *cylID, int *filID, int *cmpIDlist,
                                       int * fvecposition, int *pair_cIndex_cnIndex, double
                                       *params, int *numpairs, int *NL, int *params2) {
    const unsigned int thread_idx = (blockIdx.x * blockDim.x) + threadIdx.x;
    double cyl[3], cyln[3];
    double rmin = params[0];
    double rmax = params[1];
    int nint = params2[0];
    int fullstate = params2[1];
    bool checkstate = true;
    if (thread_idx < nint) {
        int cIndex = pair_cIndex_cnIndex[2 * thread_idx];
        int cnIndex = pair_cIndex_cnIndex[2 * thread_idx +1];
        //Get the cylinders
        for (auto i = 0; i < 3; i++) {
            cyl[i] = coord_com[3 * cIndex + i];
            cyln[i] = coord_com[3 * cnIndex + i];
        }
//        printf("%d %f %f %f %f %f %f\n", thread_idx, coord_com[3 * cIndex],coord_com[3 * cIndex + 1],coord_com[3 *
//                                                                                                                    cIndex + 2],
//               coord_com[3 * cnIndex],coord_com[3 * cnIndex + 1],coord_com[3 * cnIndex + 2]);
        if (!(fullstate) && (cylID[cIndex] <= cylID[cnIndex])) //If half list,
            // ensure that only half list is
            // calculated. This ensures that same cylinder is not added.
            checkstate = false;
//        printf("C1 %d \n", checkstate);
        if ((filID[cIndex] == filID[cnIndex])) {//If they are from the same filament, make sure they are two
            // cylinders apart.
            if (fabsf(fvecposition[cIndex] - fvecposition[cnIndex]) <= 2)
                checkstate = false;
        }
//        printf("C2 %d \n", checkstate);
        if (checkstate) {
            //check for distance between the two cylinders
            double dist = twoPointDistancemixedID(cyl, cyln, 0, 0);
//            printf("C %d %d %f %f %f\n", cIndex, cnIndex, dist, rmin, rmax);
            if (dist >= rmin && dist <= rmax) {
                int numpair_prev = atomicAdd(&numpairs[0], 1);
//                numpair_prev++;
                NL[2 * numpair_prev] = cIndex;
                NL[2 * numpair_prev + 1] = cnIndex;
//                printf("C %d %d\n", cIndex, cnIndex);
//                printf("NL %d\n",numpair_prev);
            }
        }
    }
}

__global__ void resetbitonic(int* stage){
    stage[0] = 0;//stage
    stage[1] = 0;//substage
}

__global__ void incrementbitonic(int* stage){
    if(stage[1] == 0) {
        stage[0]++;
        stage[1] = stage[0];
    }
    else
        stage[1]--;
    printf("stage %d %d\n",stage[0], stage[1]);
}
__global__ void bitonicsort(int* array, int* stage, int* nint){
    const unsigned int thread_idx = (blockIdx.x * blockDim.x) + threadIdx.x;
    int vecsize = nint[0];
    if(thread_idx < vecsize) {
        int offset = (int)pow(2.0,stage[1]);
        int orderlength = (int)pow(2,stage[1]);
        int offset2 = orderlength * (thread_idx/orderlength);
        bool ordertype = (thread_idx % (int)pow(2,stage[0] + 1))< (int)pow(2,stage[0]);
        int pos1,pos2;
       if(stage[1]!=0){
            pos1 = offset2 + thread_idx;
            pos2 = offset2 + thread_idx + offset;
        }
        else{
            pos1 = 2*thread_idx;
            pos2 = 2*thread_idx + offset;
        }
            printf("thread_idx %d ordertype %d pos %d %d orderlength %d offset2 %d\n",
                   thread_idx, ordertype, pos1, pos2, orderlength, offset2);
        if (ordertype) { //ascending
            if(array[pos1] > array[pos2]){
                int temp = array[pos1];
                array[pos1] = array[pos2];
                array[pos2] = temp;
            }
        }
        else { //descending
            if(array[pos1] < array[pos2]){
                int temp = array[pos1];
                array[pos1] = array[pos2];
                array[pos2] = temp;
            }
        }
    }
    __syncthreads();
}

//__global__ void CylinderCylinderNLCUDA(double *coord_com, int *beadSet, int *cylID, int *filID, int *cmpIDlist,
//                                       int * fvecposition, int *cylvecpospercmp, int *pair_cIndex_cmp, int
//                                        *params, int *numpairs, int *NL, int *offsetvec){
//    const unsigned int thread_idx = (blockIdx.x * blockDim.x) + threadIdx.x;
//    int cIndex = pair_cIndex_cmp[ 2 * blockIdx.x];
//    int cmpID = pair_cIndex_cmp[ 2 * blockIdx.x + 1];
//    int cnIndexstart = cylvecpospercmp[2 * cmpID];
//    int cnIndexend = cylvecpospercmp[2 * cmpID + 1];
//    double cyl[3], cyln[3];
//    //Get the cylinder
//    for(auto i = 0;i <3; i++){
//        cyl[i] = coord_com[3 * cIndex + i];
//    }
//    int offset = offsetvec[blockIdx.x];
//    int rmin = params[0];
//    int rmax = params[1];
//    int fullstate = params[2];
//    bool checkstate = false;
//    //neighboring cylinder position in vector
//    int cnIndex = threadIdx.x + cnIndexstart + offset;
//        if(cnIndex < cnIndexend && cmpIDlist[cnIndex] == cmpID){ //make sure you are adding neighbors from the
//            // compartment as mentioned in input
//            if(!(fullstate) && (cylID[cIndex] > cylID[cnIndex]) ) //If half list,
//                // ensure that only half list is
//                // calculated. This ensures that same cylinder is not added.
//                checkstate = true;
//            if((filID[cIndex] == filID[cnIndex])){//If they are from the same filament, make sure they are two
//                // cylinders apart.
//                if(fvecposition[cIndex]-fvecposition[cnIndex] <= 2)
//                    checkstate = false;
//            }
//            if(checkstate) {
//                //get coordinate
//                for (auto i = 0; i < 3; i++) {
//                    cyln[i] = coord_com[3 * cnIndex + i];
//                }
//                //check for distance between the two cylinders
//                double dist = twoPointDistancemixedID(cyl, cyln, 0, 0);
//                if(dist > rmin && dist < rmax){
//                    int numpair_prev = atomicAdd(&numpairs[0], 1);
//                    numpair_prev++;
//                    NL[2 * numpair_prev] = cIndex;
//                    NL[2 * numpair_prev + 1] = cnIndex;
//                }
//
//            }
//        }
//}
#endif
#endif //CUDA_VEC_NEIGHBORLISTIMPLCUDA_H
