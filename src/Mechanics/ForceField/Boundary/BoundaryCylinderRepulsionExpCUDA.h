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

__global__ void BoundaryCylinderRepulsionadd(double *force, double *forcecopy, int *nint) {
    const unsigned int thread_idx = (blockIdx.x * blockDim.x) + threadIdx.x;

    if (thread_idx < nint[0]) {
        for (auto i = 0; i < 3; i++) {
            force[3 * thread_idx + i] =  force[3 * thread_idx + i] + forcecopy[3 * thread_idx + i];
        }
    }
}
#endif //CUDA_VEC_BOUNDARYCYLINDERREPULSIONCUDA_H
