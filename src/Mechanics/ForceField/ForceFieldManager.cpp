
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
#ifdef CUDAACCL
    if(streamF == NULL || !(CUDAcommon::getCUDAvars().conservestreams))
        CUDAcommon::handleerror(cudaStreamCreate(&streamF));
    int nint[1]; nint[0]=CGMethod::N/3;
    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_nint, sizeof(int)));
    CUDAcommon::handleerror(cudaMemcpy(gpu_nint, nint, sizeof(int), cudaMemcpyHostToDevice));
    int THREADSPERBLOCK;
    cudaDeviceProp prop;
    cudaGetDeviceProperties(&prop, 0);
    THREADSPERBLOCK = prop.maxThreadsPerBlock;
    blocksnthreads.push_back(CGMethod::N/(3*THREADSPERBLOCK) + 1);
    if(blocksnthreads[0]==1) blocksnthreads.push_back(CGMethod::N/3);
    else blocksnthreads.push_back(THREADSPERBLOCK);
#endif
}

void ForceFieldManager::cleanupAllForceFields() {

    for(auto &ff : _forceFields)
        ff->cleanup();
#ifdef CUDAACCL
    if(!(CUDAcommon::getCUDAvars().conservestreams))
        CUDAcommon::handleerror(cudaStreamDestroy(streamF));
    if(CGMethod::N/3 > 0) {
        CUDAcommon::handleerror(cudaFree(gpu_nint));
        //Memory alloted
        //@{
//        size_t allocmem = 0;
//        allocmem += sizeof(double);
//        auto c = CUDAcommon::getCUDAvars();
//        c.memincuda -= allocmem;
//        CUDAcommon::cudavars = c;
//        std::cout<<"Total allocated memory "<<c.memincuda/1024<<endl;
//        std::cout<<"Memory allocated 0 . Memory freed "<<allocmem/1024<<endl;
        //@}
        blocksnthreads.clear();
    }
#endif
}

double ForceFieldManager::computeEnergy(double *coord, double *f, double d, bool verbose) {

    double energy = 0.0;
#ifdef CUDAACCL
    auto gU_tot = CUDAcommon::getCUDAvars().gpu_energy;
    setenergytozero<<<1,1,0,streamF>>>(gU_tot);
    CUDAcommon::handleerror(cudaStreamSynchronize(streamF));
#endif
//    std::cout<<"-------"<<endl;
    for(auto &ff : _forceFields) {
//        std::cout<<ff->getName()<<endl;
//        std::cout<<"ForceField "<<ff->getName()<<endl;
        auto tempEnergy = ff->computeEnergy(coord, f, d);
#ifdef SERIAL_CUDACROSSCHECK
        for(auto strm:CUDAcommon::getCUDAvars().streamvec)
            CUDAcommon::handleerror(cudaStreamSynchronize(*strm),"FFM","FFM.cu");
        CUDAcommon::handleerror(cudaDeviceSynchronize(),"ForceField", "ForceField");
        double cuda_energy[1];
        CUDAcommon::handleerror(cudaMemcpy(cuda_energy, CUDAcommon::cudavars.gpu_energy,  sizeof(double),
                                           cudaMemcpyDeviceToHost));
//        std::cout<<ff->getName()<<endl;
//        std::cout<<"Serial Energy "<<energy+tempEnergy<<" Cuda Energy "<<cuda_energy[0]<<endl;
#endif
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
                std::cout<<"ForceField "<<ff->getName()<<" energy "<<tempEnergy<<" d "
                        <<d<<endl;

    }
//    std::cout<<"-------"<<endl;
    return energy;
}
void ForceFieldManager::computeForces(double *coord, double *f) {
    //reset to zero

    for (int i = 0; i < CGMethod::N; i++)
        f[i] = 0.0;

#ifdef CUDAACCL
    CUDAvars cvars=CUDAcommon::getCUDAvars();
//    if(streamF == NULL || !(CUDAcommon::getCUDAvars().conservestreams))
//        CUDAcommon::handleerror(cudaStreamCreate( &streamF));
    if(cross_checkclass::Aux)
        resetForcesCUDA<<<blocksnthreads[0],blocksnthreads[1],0,streamF>>>(cvars.gpu_forceAux, gpu_nint);
    else
        resetForcesCUDA<<<blocksnthreads[0],blocksnthreads[1],0,streamF>>>(cvars.gpu_force, gpu_nint);
    CUDAcommon::handleerror(cudaStreamSynchronize(streamF));
//    if(!(CUDAcommon::getCUDAvars().conservestreams))
//        CUDAcommon::handleerror(cudaStreamDestroy(streamF));

    CUDAcommon::handleerror( cudaGetLastError() ,"resetForcesCUDA", "ForceFieldManager.cu");
#endif
    //recompute
//    double *F_i = new double[CGMethod::N];
    for(auto &ff : _forceFields) {
//                std::cout<<"ForceField "<<ff->getName()<<endl;
        ff->computeForces(coord, f);

//        CUDAcommon::handleerror(cudaDeviceSynchronize());

//        if(cross_checkclass::Aux)
//            CUDAcommon::handleerror(
//                cudaMemcpy(F_i, CUDAcommon::getCUDAvars().gpu_forceAux, 3 * Bead::getBeads().size() * sizeof
//                                   (double),
//                           cudaMemcpyDeviceToHost));
//        else
//            CUDAcommon::handleerror(
//                    cudaMemcpy(F_i, CUDAcommon::getCUDAvars().gpu_force, 3 * Bead::getBeads().size() * sizeof
//                                       (double),
//                               cudaMemcpyDeviceToHost));
//        double fmax = 0.0;
//        int id=0;
//        for (auto iter = 0; iter < Bead::getBeads().size(); iter++) {
//            if(abs(F_i[3 *iter])> fmax) {fmax = abs(F_i[3*iter]);id = iter;}
//            if(abs(F_i[3 *iter +1])> fmax) {fmax = abs(F_i[3*iter +1]);id = iter;}
//            if(abs(F_i[3 *iter +2])> fmax) {fmax = abs(F_i[3*iter +2]);id = iter;}
////            std::cout << F_i[3 * iter] << " " << F_i[3 * iter + 1] << " " << F_i[3 * iter + 2] << endl;
//        }
//        std::cout <<"Fmax "<< id<<" "<<fmax<<" "<<F_i[3 * id] << " " << F_i[3 * id + 1] << " " << F_i[3 * id + 2] <<
//                                                                                                                 endl;
    }
//    delete F_i;
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
    
}

void ForceFieldManager::resetForces() {
    
    for(auto b: Bead::getBeads()) {
        b->force.assign (3, 0); //Set force to zero;
        //std::memset((void*)(&b->loadForcesP[0]), 0, sizeof(b->loadForcesP));  //Set load force to zero;
        //std::memset((void*)(&b->loadForcesM[0]), 0, sizeof(b->loadForcesM));  //Set load force to zero;
        std::fill(b->loadForcesP.begin(), b->loadForcesP.end(), 0.0); //Set load force to zero
        std::fill(b->loadForcesM.begin(), b->loadForcesM.end(), 0.0); //Set load force to zero
    }
}

#ifdef CUDAACCL
void ForceFieldManager::CUDAcopyForces(cudaStream_t stream, double *fprev, double *f) {


//    CUDAcommon::handleerror(cudaFree(CUDAcommon::getCUDAvars().gpu_forceAux));
//    double* gpu_forceAux;
//    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_forceAux, CGMethod::N * sizeof(double)));
//    CUDAvars cvars=CUDAcommon::getCUDAvars();
//    cvars.gpu_forceAux=gpu_forceAux;
//    CUDAcommon::cudavars=cvars;

//    std::cout<<"Copyforces Number of Blocks: "<<blocksnthreads[0]<<endl;
//    std::cout<<"Threads per block: "<<blocksnthreads[1]<<endl;
    copyForcesCUDA<<<blocksnthreads[0],blocksnthreads[1],0,stream>>>(f, fprev, gpu_nint);
    CUDAcommon::handleerror( cudaGetLastError(),"copyForcesCUDA", "ForceFieldManager.cu");
}

void ForceFieldManager::assignallforcemags() {

    for(auto &ff : _forceFields)
        ff->assignforcemags();
}
#endif
