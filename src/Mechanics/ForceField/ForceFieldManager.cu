
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

#include "ForceFieldManager.h"
#include "ForceFieldManagerCUDA.h"

#include "CGMethod.h"
#include "cross_check.h"

void ForceFieldManager::vectorizeAllForceFields() {

    for(auto &ff : _forceFields)
        ff->vectorize();
}

void ForceFieldManager::cleanupAllForceFields() {

    for(auto &ff : _forceFields)
        ff->cleanup();
}

double ForceFieldManager::computeEnergy(double *coord, double *f, double d, bool verbose) {

    double energy = 0;
    for(auto &ff : _forceFields) {

        auto tempEnergy = ff->computeEnergy(coord, f, d);


        if(verbose) cout << ff->getName() << " energy = " << tempEnergy << endl;

        //if energy is infinity, exit with infinity.
        if(tempEnergy <= -1) {

            //if this is the current energy, exit ungracefully
            if(d == 0.0) {

                cout << "Energy = " << tempEnergy << endl;

                cout << "Energy of system became infinite. Try adjusting minimization parameters." << endl;
                cout << "The culprit was ... " << ff->getName() << endl;

                //get the culprit in output
                ff->whoIsCulprit();

                exit(EXIT_FAILURE);
            }
                //if this is a minimization try, just return infinity
            else return numeric_limits<double>::infinity();
        }
        else energy += tempEnergy;

    }
    return energy;
}

#ifdef CROSSCHECK
void ForceFieldManager::resetForces() {

    for(auto b: Bead::getBeads()) {
        if(cross_checkclass::Aux)
        b->forceAux.assign (3,0);
        else
        b->force.assign (3, 0); //Set force to zero;
        std::memset((void*)(&b->loadForcesP[0]), 0, sizeof(b->loadForcesP));  //Set load force to zero;
        std::memset((void*)(&b->loadForcesM[0]), 0, sizeof(b->loadForcesM));  //Set load force to zero;
    }
}
#endif
void ForceFieldManager::computeForces(double *coord, double *f) {
#ifdef CROSSCHECK
    resetForces();
#endif
    //TODO change so that you don't have to copy a vector every time during minimization.
    //reset to zero
    for (int i = 0; i < CGMethod::N; i++)
        f[i] = 0.0;

#ifdef CUDAACCL
    double* gpu_force;
    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_force, CGMethod::N * sizeof(double)));
    cudaMemcpy(gpu_force, f, CGMethod::N * sizeof(double), cudaMemcpyHostToDevice);
    CUDAvars cvars=CUDAcommon::getCUDAvars();
    if(cross_checkclass::Aux) {
        cudaFree(cvars.gpu_forceAux);
        cvars.gpu_forceAux = gpu_force;
    }
    else{
        cudaFree(cvars.gpu_force);
        cvars.gpu_force = gpu_force;
    }
    CUDAcommon::cudavars=cvars;
//    double F_i[CGMethod::N];
//    cudaMemcpy(F_i, gpu_force, CGMethod::N * sizeof(double), cudaMemcpyDeviceToHost);
//    for(auto i=0;i<CGMethod::N;i++)
//        std::cout<<F_i[i]<<" ";
//    std::cout<<endl;
#endif
    //recompute
    for(auto &ff : _forceFields) ff->computeForces(coord, f);

//#ifdef CUDAACCL
//    //TODO Remove later
//    double* gpu_force;
//    cudaMalloc((void **) &gpu_force, CGMethod::N * sizeof(double));
//    cudaMemcpy(gpu_force, f, CGMethod::N * sizeof(double), cudaMemcpyHostToDevice);
//    CUDAvars cvars=CUDAcommon::getCUDAvars();
//    cvars.gpu_force=gpu_force;
//    CUDAcommon::cudavars=cvars;
//#endif
    //WILL HAVE TO COPY AUXS AFTER THIS CALL
}

void ForceFieldManager::computeLoadForces() {

    //reset
    for (auto b: Bead::getBeads()) {

        b->loadForcesP.clear();
        b->loadForcesM.clear();
    }

    for(auto &f : _forceFields)
        f->computeLoadForces();

    //reset lfi as well
    for(auto b: Bead::getBeads()) {
        b->lfip = 0;
        b->lfim = 0;
    }
}

void ForceFieldManager::copyForces(double *fprev, double *f) {

    for (int i = 0; i < CGMethod::N; i++)
        fprev[i] = f[i];
#ifdef CUDAACCL
    //TODO Change so that the pointers to forceAux and force are exchanged and pointer to force is flushed.

    cudaFree(CUDAcommon::getCUDAvars().gpu_forceAux);
    double* gpu_forceAux;
    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_forceAux, CGMethod::N * sizeof(double)));
    CUDAvars cvars=CUDAcommon::getCUDAvars();
    cvars.gpu_forceAux=gpu_forceAux;
    CUDAcommon::cudavars=cvars;

    vector<int> blocksnthreads;
    int *gpu_nint; int nint[1]; nint[0]=CGMethod::N/3;
    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_nint, sizeof(int)));
    cudaMemcpy(gpu_nint, nint, sizeof(int), cudaMemcpyHostToDevice);
    blocksnthreads.push_back(CGMethod::N/(3*THREADSPERBLOCK) + 1);
    if(blocksnthreads[0]==1) blocksnthreads.push_back(CGMethod::N/3);
    else blocksnthreads.push_back(THREADSPERBLOCK);
    std::cout<<"Copyforces Number of Blocks: "<<blocksnthreads[0]<<endl;
    std::cout<<"Threads per block: "<<blocksnthreads[1]<<endl;
    copyForcesCUDA<<<blocksnthreads[0],blocksnthreads[1]>>>(cvars.gpu_force, cvars.gpu_forceAux, gpu_nint);
    cudaFree(gpu_nint);
#endif
}
