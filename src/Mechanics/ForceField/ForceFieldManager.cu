
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
#include <algorithm>

void ForceFieldManager::vectorizeAllForceFields() {
#ifdef CUDATIMETRACK
    chrono::high_resolution_clock::time_point tbegin, tend;
    CUDAcommon::cudatime.TvectorizeFF = 0.0;
    CUDAcommon::cudatime.TvecvectorizeFF.clear();
#endif
#ifdef CUDAACCL
    // PT1 Generate single vector of energies from all FF and add them together.
    //@{
    CUDAcommon::cudavars.offset_E=0.0;
    //@}
#endif

    for (auto &ff : _forceFields)
        ff->vectorize();

#ifdef CUDATIMETRACK
    tbegin = chrono::high_resolution_clock::now();
#endif
#ifdef CUDAACCL
    //reset offset
    if (streamF == NULL || !(CUDAcommon::getCUDAvars().conservestreams))
        CUDAcommon::handleerror(cudaStreamCreate(&streamF));
    int nint[1];
    nint[0] = CGMethod::N / 3;
    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_nint, sizeof(int)));
    CUDAcommon::handleerror(cudaMemcpyAsync(gpu_nint, nint, sizeof(int),
                                        cudaMemcpyHostToDevice, streamF));
    CUDAcommon::handleerror(cudaMalloc((void **) &(CUDAcommon::cudavars.gpu_energyvec),
                                       CUDAcommon::cudavars.offset_E * sizeof(floatingpoint)));
    int THREADSPERBLOCK;
    cudaDeviceProp prop;
    cudaGetDeviceProperties(&prop, 0);
    THREADSPERBLOCK = prop.maxThreadsPerBlock;

    blocksnthreads.push_back(CGMethod::N / (3 * THREADSPERBLOCK) + 1);
    if (blocksnthreads[0] == 1) blocksnthreads.push_back(CGMethod::N / 3);
    else blocksnthreads.push_back(THREADSPERBLOCK);

    // PT2 Generate single vector of energies from all FF and add them together.
    //@{
//    std::cout<<"CUDA energy total nint "<<CUDAcommon::cudavars.offset_E<<endl;
    bntaddvec2.clear();
    bntaddvec2 = getaddred2bnt(CUDAcommon::cudavars.offset_E);
    CUDAcommon::handleerror(cudaMalloc((void **) &(CUDAcommon::cudavars.gpu_energyvec), bntaddvec2.at
            (0)*sizeof (floatingpoint)));
    vector<floatingpoint> zerovec(bntaddvec2.at(0));
    fill(zerovec.begin(),zerovec.begin()+bntaddvec2.at(0),0.0);
    CUDAcommon::handleerror(cudaMemcpyAsync(CUDAcommon::cudavars.gpu_energyvec, zerovec.data(),
                            bntaddvec2.at(0) * sizeof(floatingpoint), cudaMemcpyHostToDevice,streamF));
/*    CUDAcommon::handleerror(cudaMemsetAsync(CUDAcommon::cudavars.gpu_energyvec, 0,
                                            bntaddvec2.at(0) * sizeof(floatingpoint), streamF));*/

    params.clear();
    params.push_back(CUDAcommon::cudavars.offset_E);
    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_params, sizeof(int)));
    CUDAcommon::handleerror(cudaMemcpyAsync(gpu_params, params.data(), sizeof(int),
                                       cudaMemcpyHostToDevice, streamF));
    //@}
#endif
#ifdef CUDATIMETRACK
//    CUDAcommon::handleerror(cudaDeviceSynchronize(),"CGPolakRibiereMethod.cu","CGPolakRibiereMethod.cu");
    tend= chrono::high_resolution_clock::now();
    chrono::duration<floatingpoint> elapsed_run(tend - tbegin);
    CUDAcommon::cudatime.TvectorizeFF += elapsed_run.count();
    std::cout<<"Time total vectorizeFF (s) "<<CUDAcommon::cudatime.TvectorizeFF<<endl;
    std::cout<<"Time split vectorizeFF (s) ";
    for(auto x:CUDAcommon::cudatime.TvecvectorizeFF)
        std::cout<<x<<" ";
    std::cout<<endl;
#endif
}

void ForceFieldManager::cleanupAllForceFields() {

    for (auto &ff : _forceFields)
        ff->cleanup();
#ifdef CUDAACCL
    if (!(CUDAcommon::getCUDAvars().conservestreams))
        CUDAcommon::handleerror(cudaStreamDestroy(streamF));
    //cleanup energy vector
    // single vector of energies from all FF and add them together.
    //@{
    CUDAcommon::handleerror(cudaFree(CUDAcommon::cudavars.gpu_energyvec));
    CUDAcommon::handleerror(cudaFree(gpu_params));
    //@}
    if (CGMethod::N / 3 > 0) {
        CUDAcommon::handleerror(cudaFree(gpu_nint));
        //Memory alloted
        //@{
//        size_t allocmem = 0;
//        allocmem += sizeof(floatingpoint);
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

floatingpoint ForceFieldManager::computeEnergy(floatingpoint *coord, floatingpoint *f, floatingpoint d, bool verbose) {
#ifdef CUDATIMETRACK
    chrono::high_resolution_clock::time_point tbegin, tend;
//    CUDAcommon::cudatime.TcomputeE = 0.0;
    CUDAcommon::cudatime.TveccomputeE.clear();
    CUDAcommon::cudatime.Ecount++;
//    CUDAcommon::serltime.TcomputeE = 0.0;
    CUDAcommon::serltime.TveccomputeE.clear();
    CUDAcommon::serltime.Ecount++;
#endif
    floatingpoint energy = 0.0;
#ifdef CUDAACCL
#ifdef CUDA_INDIVIDUAL_ESUM
    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_Uvec, sizeof (floatingpoint)));
    CUDAcommon::handleerror(cudaMemset(gpu_Uvec, 0.0, sizeof (floatingpoint)));
#else
    floatingpoint *gpu_Uvec = CUDAcommon::getCUDAvars().gpu_energy;
    /*floatingpoint *gpu_Uvec;
    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_Uvec, sizeof (floatingpoint)));
    CUDAcommon::handleerror(cudaMemsetAsync(gpu_Uvec, 0, sizeof (floatingpoint),streamF));*/
#endif
#ifdef SERIAL_CUDACROSSCHECK
    CUDAcommon::handleerror(cudaMemset(CUDAcommon::cudavars.gpu_energyvec, 0, bntaddvec2.at(0) * sizeof
            (floatingpoint)));
#endif
#ifdef CUDATIMETRACK
    tbegin = chrono::high_resolution_clock::now();
#endif
/*    auto gU_tot = CUDAcommon::getCUDAvars().gpu_energy;
    setenergytozero << < 1, 1, 0, streamF >> > (gU_tot);*/
    CUDAcommon::handleerror(cudaStreamSynchronize(streamF));
#ifdef CUDATIMETRACK
//    CUDAcommon::handleerror(cudaDeviceSynchronize(),"ForceFieldManager.cu",
//                            "computeEnergy");
    tend= chrono::high_resolution_clock::now();
    chrono::duration<floatingpoint> elapsed_run(tend - tbegin);
    CUDAcommon::cudatime.TcomputeE += elapsed_run.count();
    CUDAcommon::cudatime.TcomputeEiter += elapsed_run.count();
#endif
#endif
#ifdef SERIAL_CUDACROSSCHECK
    CUDAcommon::handleerror(cudaDeviceSynchronize());
    floatingpoint cuda_lambda[1];
    CUDAcommon::handleerror(cudaMemcpy(cuda_lambda, CUDAcommon::cudavars.gpu_lambda,  sizeof(floatingpoint),
                                       cudaMemcpyDeviceToHost));

    std::cout<<"Lambda used CUDA "<<cuda_lambda[0]<<" SERL "<<d<<endl;
#endif
    for (auto &ff : _forceFields) {

        auto tempEnergy = ff->computeEnergy(coord, f, d);
#ifdef ALLSYNC
        cudaDeviceSynchronize();
#endif
//        std::cout<<ff->getName()<<" "<<tempEnergy<<endl;
        if (verbose) cout << ff->getName() << " energy = " << tempEnergy << endl;
        //if energy is infinity, exit with infinity.
        if (tempEnergy <= -1) {

            //if this is the current energy, exit ungracefully
            if (d == 0.0) {

                cout << "Energy = " << tempEnergy << endl;

                cout
                        << "Energy of system became infinite. Try adjusting minimization parameters."
                        << endl;
                cout << "The culprit was ... " << ff->getName() << endl;

                //get the culprit in output
                ff->whoIsCulprit();

                exit(EXIT_FAILURE);
            }
                //if this is a minimization try, just return infinity
            else return numeric_limits<floatingpoint>::infinity();
        } else energy += tempEnergy;
#ifdef SERIAL_CUDACROSSCHECK
        cudaDeviceSynchronize();
        resetfloatingpointvariableCUDA<<<1,1,0, streamF>>>(gpu_Uvec);
        addvectorred3<<<bntaddvec2.at(2),bntaddvec2.at(3), bntaddvec2.at(3) * sizeof(floatingpoint)
                                                                   , streamF>>>
        (CUDAcommon::cudavars.gpu_energyvec, gpu_params,
                gpu_Uvec);
        floatingpoint cuda_energyvec[1];
        CUDAcommon::handleerror(cudaMemcpy(cuda_energyvec, gpu_Uvec, sizeof(floatingpoint),
                                           cudaMemcpyDeviceToHost));
        std::cout<<ff->getName()<<" Energy. CUDA "<<cuda_energyvec[0]<<" SERL "
                ""<<energy<<endl;
#endif
    }
#ifdef CUDATIMETRACK
    tbegin = chrono::high_resolution_clock::now();
#endif
    //Add energies
#ifdef CUDAACCL
//    std::cout<<"Total nint "<<bntaddvec2.at(0)<<" "<<CUDAcommon::cudavars.offset_E<<endl;
    //Synchronize streams
    for(auto strm:CUDAcommon::getCUDAvars().streamvec) {
            CUDAcommon::handleerror(cudaStreamSynchronize(*strm), "computeEnergy",
                                    "ForceFieldManager.cu");
        }
#ifdef CUDATIMETRACK
    tend= chrono::high_resolution_clock::now();
    chrono::duration<floatingpoint> elapsed_run2(tend - tbegin);
    CUDAcommon::cudatime.TveccomputeE.push_back(elapsed_run2.count());
    CUDAcommon::cudatime.TcomputeE += elapsed_run2.count();
    CUDAcommon::cudatime.TcomputeEiter += elapsed_run2.count();
    tbegin = chrono::high_resolution_clock::now();
#endif
//    std::cout<<"CUDA energy total nint "<<CUDAcommon::cudavars.offset_E<<endl;
    /*vector<floatingpoint> ones;
    for(int i = 0;i<8192;i++)
        ones.push_back(1.0);
    CUDAcommon::handleerror(cudaMemcpyAsync(CUDAcommon::cudavars.gpu_energyvec, ones
                                                        .data() ,
                                                bntaddvec2.at(0) * sizeof
                                                        (floatingpoint),
                                                cudaMemcpyHostToDevice,streamF));*/
/*    CUDAcommon::handleerror(cudaMemsetAsync(CUDAcommon::cudavars.gpu_energyvec, 1,
                                            bntaddvec2.at(0) * sizeof
            (floatingpoint),streamF));
    cudaDeviceSynchronize();*/
    resetfloatingpointvariableCUDA<<<1,1,0, streamF>>>(gpu_Uvec);
    addvectorred3<<<bntaddvec2.at(2),bntaddvec2.at(3), bntaddvec2.at(3) * sizeof(floatingpoint)
                    , streamF>>>(CUDAcommon::cudavars.gpu_energyvec, gpu_params,
                        gpu_Uvec);
    CUDAcommon::handleerror(cudaStreamSynchronize(streamF));
#ifdef DETAILEDOUTPUT_ENERGY
    cudaDeviceSynchronize();
    floatingpoint cuda_energyvec[1];
    CUDAcommon::handleerror(cudaMemcpy(cuda_energyvec, gpu_Uvec, sizeof(floatingpoint),
                                       cudaMemcpyDeviceToHost));
    std::cout<<"vector energy addition CUDA "<<cuda_energyvec[0]<<" SERL "<<energy<<endl;
#endif
//    CUDAcommon::handleerror(cudaFree(CUDAcommon::cudavars.gpu_energyvec));
#endif
#ifdef ALLSYNC
    cudaDeviceSynchronize();
#endif

#ifdef CUDATIMETRACK
    tend= chrono::high_resolution_clock::now();
    chrono::duration<floatingpoint> elapsed_run3(tend - tbegin);
    CUDAcommon::cudatime.TveccomputeE.push_back(elapsed_run3.count());
    CUDAcommon::cudatime.TcomputeE += elapsed_run3.count();
    CUDAcommon::cudatime.TcomputeEiter += elapsed_run3.count();
//    std::cout<<"Time total computeEnergy (s) CUDA "<<CUDAcommon::cudatime
//            .TcomputeE<<" SERL "<<CUDAcommon::serltime.TcomputeE<<" factor "
//                     ""<<CUDAcommon::serltime.TcomputeE/CUDAcommon::cudatime.TcomputeE<<endl;
//    std::cout<<"Time split computeEnergy (s) CUDA ";
//    for(auto x:CUDAcommon::cudatime.TveccomputeE)
//        std::cout<<x<<" ";
//    std::cout<<endl;
//    std::cout<<"Time split computeEnergy (s) SERL ";
//    for(auto x:CUDAcommon::serltime.TveccomputeE)
//        std::cout<<x<<" ";
//    std::cout<<endl;
#endif
    return energy;
}

void ForceFieldManager::computeForces(floatingpoint *coord, floatingpoint *f) {
    //reset to zero
#ifdef CUDATIMETRACK
    chrono::high_resolution_clock::time_point tbegin, tend;
    CUDAcommon::cudatime.TcomputeF = 0.0;
    CUDAcommon::cudatime.TveccomputeF.clear();
    CUDAcommon::serltime.TcomputeF = 0.0;
    CUDAcommon::serltime.TveccomputeF.clear();
    tbegin = chrono::high_resolution_clock::now();
#endif
    //@{
    for (int i = 0; i < CGMethod::N; i++)
        f[i] = 0.0;
    //@}
#ifdef CUDATIMETRACK
    tend= chrono::high_resolution_clock::now();
    chrono::duration<floatingpoint> elapsed_run(tend - tbegin);
    CUDAcommon::serltime.TveccomputeF.push_back(elapsed_run.count());
    CUDAcommon::serltime.TcomputeF += elapsed_run.count();
    tbegin = chrono::high_resolution_clock::now();
#endif
#ifdef CUDAACCL
    CUDAvars cvars = CUDAcommon::getCUDAvars();
    if (cross_checkclass::Aux)
        resetForcesCUDA << < blocksnthreads[0], blocksnthreads[1], 0, streamF >> >
                                                                      (cvars.gpu_forceAux, gpu_nint);
    else
        resetForcesCUDA << < blocksnthreads[0], blocksnthreads[1], 0, streamF >> >
                                                                      (cvars.gpu_force, gpu_nint);
    CUDAcommon::handleerror(cudaStreamSynchronize(streamF));

    CUDAcommon::handleerror(cudaGetLastError(), "resetForcesCUDA", "ForceFieldManager.cu");
#endif
#ifdef CUDATIMETRACK
    tend= chrono::high_resolution_clock::now();
    chrono::duration<floatingpoint> elapsed_run2(tend - tbegin);
    CUDAcommon::cudatime.TveccomputeF.push_back(elapsed_run2.count());
    CUDAcommon::cudatime.TcomputeF += elapsed_run2.count();
    tbegin = chrono::high_resolution_clock::now();
#endif
    //recompute
//    floatingpoint *F_i = new floatingpoint[CGMethod::N];
    for (auto &ff : _forceFields) {
        ff->computeForces(coord, f);
#ifdef ALLSYNC
        cudaDeviceSynchronize();
#endif

//        if(cross_checkclass::Aux)
//            CUDAcommon::handleerror(
//                cudaMemcpy(F_i, CUDAcommon::getCUDAvars().gpu_forceAux, 3 * Bead::getBeads().size() * sizeof
//                                   (floatingpoint),
//                           cudaMemcpyDeviceToHost));
//        else
//            CUDAcommon::handleerror(
//                    cudaMemcpy(F_i, CUDAcommon::getCUDAvars().gpu_force, 3 * Bead::getBeads().size() * sizeof
//                                       (floatingpoint),
//                               cudaMemcpyDeviceToHost));
//        floatingpoint fmax = 0.0;
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
        std::fill(b->loadForcesM.begin(), b->loadForcesM.end(), 0.0);
        std::fill(b->loadForcesP.begin(), b->loadForcesP.end(), 0.0);
//        b->loadForcesP.clear();
//        b->loadForcesM.clear();
    }

    for (auto &f : _forceFields)
        f->computeLoadForces();

    //reset lfi as well
    for (auto b: Bead::getBeads()) {
        b->lfip = 0;
        b->lfim = 0;
    }
}

void ForceFieldManager::copyForces(floatingpoint *fprev, floatingpoint *f) {

    for (int i = 0; i < CGMethod::N; i++)
        fprev[i] = f[i];
}

#ifdef CUDAACCL

void ForceFieldManager::CUDAcopyForces(cudaStream_t stream, floatingpoint *fprev, floatingpoint *f) {


//    CUDAcommon::handleerror(cudaFree(CUDAcommon::getCUDAvars().gpu_forceAux));
//    floatingpoint* gpu_forceAux;
//    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_forceAux, CGMethod::N * sizeof(floatingpoint)));
//    CUDAvars cvars=CUDAcommon::getCUDAvars();
//    cvars.gpu_forceAux=gpu_forceAux;
//    CUDAcommon::cudavars=cvars;

//    std::cout<<"Copyforces Number of Blocks: "<<blocksnthreads[0]<<endl;
//    std::cout<<"Threads per block: "<<blocksnthreads[1]<<endl;
    copyForcesCUDA << < blocksnthreads[0], blocksnthreads[1], 0, stream >> >
                                                                 (f, fprev, gpu_nint);
    CUDAcommon::handleerror(cudaGetLastError(), "copyForcesCUDA", "ForceFieldManager.cu");
}

void ForceFieldManager::assignallforcemags() {

    for (auto &ff : _forceFields)
        ff->assignforcemags();
}

#endif