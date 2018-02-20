//
// Created by aravind on 11/1/17.
//

#ifndef CUDA_VEC_CGMETHODCUDA_H
#define CUDA_VEC_CGMETHODCUDA_H
#include <cuda.h>
#include <cuda_runtime.h>
#include "math.h"
#include "utility.h"
#include "assert.h"

__global__ void moveBeadsCUDA(double *coord, double* f, double *d,  int *nint, bool *checkin) {
    if(checkin[0] == false) return; //if it is not in minimization state
    const unsigned int thread_idx = (blockIdx.x * blockDim.x) + threadIdx.x;
//    if(thread_idx == 0)
//        printf("%d %f\n", checkin[0], d[0]);

    if (thread_idx < nint[0]) {
        for (auto i = 0; i < 3; i++) {
            coord[3 * thread_idx + i] = coord[3 * thread_idx + i] + d[0] * f[3 * thread_idx + i] ;
        }
    }
}
__global__ void shiftGradientCUDA(double *f, double* fAux, int * nint, double* newGrad, double* prevGrad, double*
curGrad, bool *Mstate) {
    if(Mstate[0] == false) return;
    const unsigned int thread_idx = (blockIdx.x * blockDim.x) + threadIdx.x;
    double d  = fmax(0.0, (newGrad[0] - prevGrad[0]) / curGrad[0]);
//    if(thread_idx == 0)
//        printf("Shift Gradient %f %f %f %f\n", d, newGrad[0], prevGrad[0],curGrad[0] );

    if (thread_idx < nint[0]) {
        for (auto i = 0; i < 3; i++) {
            f[3 * thread_idx + i] = fAux[3 * thread_idx + i] + d * f[3 * thread_idx + i];
        }
    }
}

__global__ void shiftGradientCUDAifsafe(double *f, double* fAux, int * nint, bool *Mstate, bool *Sstate) {
    if(Mstate[0] == false || Sstate[0] == false) return;//checks for Minimization state and Safe state.
    const unsigned int thread_idx = (blockIdx.x * blockDim.x) + threadIdx.x;
//    if(thread_idx == 0)
//        printf("shiftGradient safe \n");
    if (thread_idx < nint[0]) {
        for (auto i = 0; i < 3; i++) {
            f[3 * thread_idx + i] = fAux[3 * thread_idx + i];
        }
    }
}

__global__ void allFADotFCUDA(double *f1, double *f2, double *g, int * nint) {
    const unsigned int thread_idx = (blockIdx.x * blockDim.x) + threadIdx.x;
    if (thread_idx < nint[0]) {
            g[thread_idx] = 0.0;
    }
    __syncthreads();
    if (thread_idx < nint[0]) {
        for (auto i = 0; i < 3; i++) {
           g[thread_idx] +=f1[3 * thread_idx +i] * f2[3 * thread_idx +i];
        }
    }
}

__global__ void maxFCUDA(double *f, int * nint, double *fmax) {

    double mag;
    fmax[0]=0.0;
    for(int i = 0; i < 3 * nint[0]; i++) {
        mag = sqrt(f[i]*f[i]);
        if(mag > fmax[0]) {fmax[0] = mag;}
    }
//    id = id -id%3;
//    printf("Fmax %d %f %f %f %f\n", id, fmax[0], f[id], f[id+1],f[id+2]);
}

__global__ void addvector(double *U, int *params, double *U_sum){
    U_sum[0] = 0.0;
//    printf("%d \n", params[0]);
    for(auto i=0;i<params[0];i++){
        U_sum[0]  += U[i];
    }
}

__global__ void initializeLambdaCUDA(bool *checkin, bool* checkout, double *currentEnergy, double *energy,
                                     double* CUDA_lambda, double* fmax, double* params, bool *Safestate){
    checkin[0] = false;
    checkout[0] = false;
    double LAMBDAMAX = params[3];
//    printf("SS %d \n",Safestate[0]);
    if(Safestate[0] == true){//safebacktrackinglinesearch
        CUDA_lambda[0] = LAMBDAMAX;
    }
    else{//backtrackinglinesearch
        double MAXDIST = params[4];

        if(fmax[0]==0.0) {
            CUDA_lambda[0] = 0.0;
            checkout[0] = true;
        }
        else
            CUDA_lambda[0] = fmin(LAMBDAMAX, MAXDIST / fmax[0]);
    }
//    printf("CL %f %f\n", LAMBDAMAX, CUDA_lambda[0]);
}

__global__ void resetlambdaCUDA (double *CUDA_lambda){
    CUDA_lambda[0] = 0.0;
}
//__global__ void prepbacktracking(bool *checkin, bool* checkout, double *currentEnergy, double *energy,
//                                 double* CUDA_lambda, double* fmax, double* params){
//    //if(checkin[0]) return;
//    checkin[0] = false;
//    checkout[0] = false;
//    currentEnergy[0] = energy[0];
//    checkout[0] = false;
//    double LAMBDAMAX = params[3];
//    double MAXDIST = params[4];
//
//    if(fmax[0]==0.0) {
//       CUDA_lambda[0] = 0.0;
//        checkout[0] = true;
//    }
//    else
//        CUDA_lambda[0] = fmin(LAMBDAMAX, MAXDIST / fmax[0]);
//    //printf("%f \n", CUDA_lambda[0]);
//}
//
//__global__ void prepsafebacktracking(bool* checkin, bool *checkout, double* currentEnergy, double* energy, double*
//CUDA_lambda, double* params) {
//    checkin[0] = false;
//    checkout[0] = false;
//    currentEnergy[0] = energy[0];
//    double LAMBDAMAX = params[3];
//    CUDA_lambda[0] = LAMBDAMAX;
//}
__global__ void setcurrentenergy( double* energy, double* currentenergy){
    currentenergy[0] = energy[0];
}

__global__ void findLambdaCUDA(double* energyLambda, double* currentEnergy, double* FDotFA, double *fmax, double*
lambda, double *params, bool* prev_convergence, bool*  current_convergence, bool *safestate){
    if(prev_convergence[0]) return;
    current_convergence[0] = false;
    double LAMBDAREDUCE = params[1];
    double LAMBDATOL = params[2];
    double MAXDIST = params[4];
    double idealEnergyChange = 0.0;
    double energyChange = 0.0;
//    printf("Energies %f %f\n", energyLambda[0],currentEnergy[0]);
    if(safestate[0] == true){
        energyChange = energyLambda[0] - currentEnergy[0];
        if (energyChange <= 0.0) {
            current_convergence[0] = true;
//            printf("lambda_converged %.8f %.8f\n", lambda[0], LAMBDATOL);
            return;
        }
    }
    else {
        double BACKTRACKSLOPE = params[0];
        idealEnergyChange = -BACKTRACKSLOPE * lambda[0] * FDotFA[0];
        energyChange = energyLambda[0] - currentEnergy[0];

        if (energyChange <= idealEnergyChange) {
            current_convergence[0] = true;
//            printf("lambda_converged %.8f %.8f\n", lambda[0], LAMBDATOL);
            return;
        }
    }
    lambda[0] *= LAMBDAREDUCE;

    if(lambda[0] <= 0.0 || lambda[0] <= LAMBDATOL) {
        current_convergence[0] = true;
        if(safestate[0] == true)
            lambda[0] = MAXDIST / fmax[0];
        else
            lambda[0] = 0.0;
    }
//    printf("lambda %.8f %.8f\n", lambda[0], LAMBDATOL);
}
//__global__ void CUDAbacktrackingfindlambda(double* energyLambda, double* currentEnergy, double* FDotFA, double*
//lambda, double *params, bool* prev_convergence, bool*  current_convergence){
////    printf("%d %d %f %f\n", prev_convergence[0],current_convergence[0],energyLambda[0],currentEnergy[0]);
//    if(prev_convergence[0]) return;
//    current_convergence[0] = false;
//    double BACKTRACKSLOPE = params[0];
//    double LAMBDAREDUCE = params[1];
//    double LAMBDATOL = params[2];
//    double idealEnergyChange = -BACKTRACKSLOPE * lambda[0] * FDotFA[0];
//    double energyChange =  energyLambda[0] - currentEnergy[0];
////    printf("%f \n", lambda[0]);
//    if(energyChange <= idealEnergyChange) {
//        current_convergence[0] = true;
//        return;}
//    lambda[0] *= LAMBDAREDUCE;
////    printf("lambdareduce %f \n", lambda[0]);
//    if(lambda[0] <= 0.0 || lambda[0] <= LAMBDATOL) {
//        current_convergence[0] = true;
//        lambda[0] = 0.0;
//    }
//}

//__global__ void CUDAsafebacktrackingfindlambda(double* energyLambda, double* currentenergy, double* fmax,
//                                               double* lambda, double* params,  bool* prev_convergence,
//                                               bool*  current_convergence){
//    if(prev_convergence[0]) return;
//    current_convergence[0] = false;
//    double energyChange = energyLambda[0] - currentenergy[0];
//    double LAMBDAREDUCE = params[1];
//    double LAMBDATOL = params[2];
//    double MAXDIST = params[4];
//    //return if ok
//    if(energyChange <= 0.0) {current_convergence[0] = true; return;}
//    //reduce lambda
//    lambda[0] *= LAMBDAREDUCE;
//
//    //just shake if we cant find an energy min,
//    //so we dont get stuck
//    if(lambda[0] <= 0.0 || lambda[0] <= LAMBDATOL) {current_convergence[0] = true; lambda[0] = MAXDIST / fmax[0];  }
//}

#ifdef CUDAACCL

__global__ void getsafestateCUDA(double* FDotFA, double* curGrad, double* newGrad, bool* checkout){
    //if(checkin[0] == true) return;
    checkout[0] = false;
    if(FDotFA[0] <= 0.0 || areEqual(curGrad[0], newGrad[0]))
        checkout[0] = true;
//    printf("safe state %d \n", checkout[0]);
    curGrad[0] = newGrad[0];
}
__global__ void getminimizestateCUDA(double *fmax, double *GRADTOL, bool *checkin, bool *checkout) {
    //maxF() > GRADTOL
//    printf("minstate %f %f %d %d\n", fmax[0],GRADTOL[0],checkin[0],checkout[0]);
    checkout[0] = false;
    if(checkin[0] == false) return;
    if(fmax[0] > GRADTOL[0])
        checkout[0] = true;
//    printf("minstate %f %f %d %d\n", fmax[0],GRADTOL[0],checkin[0],checkout[0]);
}

__global__ void initializePolak(bool* Mcheckin, bool *Mcheckout, bool *Scheckin, bool *Scheckout){
    Mcheckin[0] = true;
    Mcheckout[0] = true;
    Scheckin[0] = false;
    Scheckout[0] = false;
}
#endif
#endif //CUDA_VEC_CGMETHODCUDA_H
