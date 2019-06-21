//
// Created by aravind on 11/21/17.
//

#ifndef CUDA_VEC_BOUNDARYCYLINDERREPULSIONCUDA_H
#define CUDA_VEC_BOUNDARYCYLINDERREPULSIONCUDA_H
#ifdef CUDAACCL
#include "MathFunctions.h"
#include "SysParams.h"
#include <limits>
#include <assert.h>
#include <cuda.h>
#include <cuda_runtime.h>
using namespace mathfunc;


//#if !defined(__CUDA_ARCH__) || __CUDA_ARCH__ >= 600
//
//#else
//static __inline__ __device__ floatingpoint atomicAdd(floatingpoint *address, floatingpoint val) {
//    unsigned long long int* address_as_ull = (unsigned long long int*)address;
//    unsigned long long int old = *address_as_ull, assumed;
//    if (val==0.0)
//      return __longlong_as_floatingpoint(old);
//    do {
//      assumed = old;
//      old = atomicCAS(address_as_ull, assumed, __floatingpoint_as_longlong(val +__longlong_as_floatingpoint(assumed)));
//    } while (assumed != old);
//    return __longlong_as_floatingpoint(old);
//  }
//
//#endif
__global__ void BoundaryCylinderRepulsionE(floatingpoint *U_b, floatingpoint *U_tot) {

    atomicAdd(&U_tot[0], U_b[0]);
}
#endif
#endif //CUDA_VEC_BOUNDARYCYLINDERREPULSIONCUDA_H
