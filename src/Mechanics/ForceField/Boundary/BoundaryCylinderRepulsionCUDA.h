//
// Created by aravind on 11/21/17.
//

#ifndef CUDA_VEC_BOUNDARYCYLINDERREPULSIONCUDA_H
#define CUDA_VEC_BOUNDARYCYLINDERREPULSIONCUDA_H

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
//#endif
__global__ void BoundaryCylinderRepulsionE(double *U_b, double *U_tot) {

    atomicAdd(&U_tot[0], U_b[0]);
}
#endif
#endif //CUDA_VEC_BOUNDARYCYLINDERREPULSIONCUDA_H
