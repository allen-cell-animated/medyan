//
// Created by aravind on 11/1/17.
//

#ifndef CUDA_VEC_CGMETHODCUDA_H
#define CUDA_VEC_CGMETHODCUDA_H
#include <cuda.h>
#include <cuda_runtime.h>

__global__ void moveBeadsgpu(double *coord, double* f, double *d,  int * nint) {
    const unsigned int thread_idx = (blockIdx.x * blockDim.x) + threadIdx.x;

    if (thread_idx < nint[0]) {
        for (auto i = 0; i < 3; i++) {
            coord[3 * thread_idx + i] = coord[3 * thread_idx + i] + d[0] * f[3 * thread_idx + i] ;
        }
    }
}

__global__ void shiftGradientCUDA(double *f, double* fAux, double *d,  int * nint) {
    const unsigned int thread_idx = (blockIdx.x * blockDim.x) + threadIdx.x;

    if (thread_idx < nint[0]) {
        for (auto i = 0; i < 3; i++) {
            f[3 * thread_idx + i] = fAux[3 * thread_idx + i] + d[0] * f[3 * thread_idx + i];
        }
    }
}
#endif //CUDA_VEC_CGMETHODCUDA_H
