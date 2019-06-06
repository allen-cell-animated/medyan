
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v3.2.1
//
//  Copyright (2015-2018)  Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the MEDYAN web page for more information:
//  http://www.medyan.org
//------------------------------------------------------------------

#include "ForceFieldManager.h"

#include <numeric> // iota

#include "ForceFieldManagerCUDA.h"

#include "CGMethod.h"
#include "cross_check.h"
#include <algorithm>

#include "SubSystem.h"
#include "Structure/Bead.h"
#include "Structure/SurfaceMesh/Membrane.hpp"
#include "VisualHelper.hpp"

namespace {

template< bool stretched > void updateGeometryValue() {
    for(auto m : Membrane::getMembranes()) m->updateGeometryValue<stretched>();
}

void updateGeometryValueWithDerivative() {
    for(auto m : Membrane::getMembranes()) m->updateGeometryValueWithDerivative();
}

void prepareForceSharedData() {
    // TODO
    // std::lock_guard<std::mutex> guard(visual::shared::dataMutex);

    // visual::shared::arrowVertexCoords.resize(2 * Bead::getDbDataConst().forces.size_raw());
    // visual::shared::lineVertexIndices.resize(2 * Bead::numBeads());
    // std::iota(visual::shared::lineVertexIndices.begin(), visual::shared::lineVertexIndices.end(), 0u);

    // visual::shared::forceChanged = true;
    // visual::shared::forceIndexChanged = true;
}

void updateForceSharedData() {
    // TODO
    // std::lock_guard<std::mutex> guard(visual::shared::dataMutex);

    // constexpr float factor = 5.0f;

    // size_t numBeads = Bead::numBeads();
    // for(size_t i = 0; i < numBeads; ++i) {
    //     const auto coord = Bead::getDbDataConst().coords[i];
    //     const auto force = Bead::getDbDataConst().forces[i];
    //     const auto endCoord = coord + factor * force;
    //     visual::shared::arrowVertexCoords[6 * i] = coord[0];
    //     visual::shared::arrowVertexCoords[6 * i + 1] = coord[1];
    //     visual::shared::arrowVertexCoords[6 * i + 2] = coord[2];
    //     visual::shared::arrowVertexCoords[6 * i + 3] = endCoord[0];
    //     visual::shared::arrowVertexCoords[6 * i + 4] = endCoord[1];
    //     visual::shared::arrowVertexCoords[6 * i + 5] = endCoord[2];
    // }

    // visual::shared::forceChanged = true;
}

template< bool stretched >
void updateMembraneSharedData() {
    // TODO stretched version
    if(!stretched) {
        visual::copySystemDataAndRunHelper(visual::sys_data_update::BeadPosition);
    }
}

} // namespace

void ForceFieldManager::vectorizeAllForceFields() {
    prepareForceSharedData();
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
                                       CUDAcommon::cudavars.offset_E * sizeof(double)));
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
            (0)*sizeof (double)));
    vector<double> zerovec(bntaddvec2.at(0));
    fill(zerovec.begin(),zerovec.begin()+bntaddvec2.at(0),0.0);
    CUDAcommon::handleerror(cudaMemcpyAsync(CUDAcommon::cudavars.gpu_energyvec, zerovec.data(),
                            bntaddvec2.at(0) * sizeof(double), cudaMemcpyHostToDevice,streamF));
/*    CUDAcommon::handleerror(cudaMemsetAsync(CUDAcommon::cudavars.gpu_energyvec, 0,
                                            bntaddvec2.at(0) * sizeof(double), streamF));*/

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
    chrono::duration<double> elapsed_run(tend - tbegin);
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

template< bool stretched >
double ForceFieldManager::computeEnergy(double *coord, bool verbose) const {
#ifdef CUDATIMETRACK
    chrono::high_resolution_clock::time_point tbegin, tend;
//    CUDAcommon::cudatime.TcomputeE = 0.0;
    CUDAcommon::cudatime.TveccomputeE.clear();
    CUDAcommon::cudatime.Ecount++;
//    CUDAcommon::serltime.TcomputeE = 0.0;
    CUDAcommon::serltime.TveccomputeE.clear();
    CUDAcommon::serltime.Ecount++;
#endif
    double energy = 0.0;
#ifdef CUDAACCL
#ifdef CUDA_INDIVIDUAL_ESUM
    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_Uvec, sizeof (double)));
    CUDAcommon::handleerror(cudaMemset(gpu_Uvec, 0.0, sizeof (double)));
#else
    double *gpu_Uvec = CUDAcommon::getCUDAvars().gpu_energy;
    /*double *gpu_Uvec;
    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_Uvec, sizeof (double)));
    CUDAcommon::handleerror(cudaMemsetAsync(gpu_Uvec, 0, sizeof (double),streamF));*/
#endif
#ifdef SERIAL_CUDACROSSCHECK
    CUDAcommon::handleerror(cudaMemset(CUDAcommon::cudavars.gpu_energyvec, 0, bntaddvec2.at(0) * sizeof
            (double)));
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
    chrono::duration<double> elapsed_run(tend - tbegin);
    CUDAcommon::cudatime.TcomputeE += elapsed_run.count();
    CUDAcommon::cudatime.TcomputeEiter += elapsed_run.count();
#endif
#endif
#ifdef SERIAL_CUDACROSSCHECK
    CUDAcommon::handleerror(cudaDeviceSynchronize());
    double cuda_lambda[1];
    CUDAcommon::handleerror(cudaMemcpy(cuda_lambda, CUDAcommon::cudavars.gpu_lambda,  sizeof(double),
                                       cudaMemcpyDeviceToHost));

    std::cout<<"Lambda used CUDA "<<cuda_lambda[0]<<" SERL "<<d<<endl;
#endif

    // Compute geometry
    updateGeometryValue< stretched >();

    for (auto &ff : _forceFields) {

        auto tempEnergy = ff->computeEnergy(coord, stretched);
#ifdef ALLSYNC
        cudaDeviceSynchronize();
#endif
//        std::cout<<ff->getName()<<" "<<tempEnergy<<endl;
        if (verbose) cout << ff->getName() << " energy = " << tempEnergy << endl;
        //if energy is infinity, exit with infinity.
        if (tempEnergy <= -1) {

            //if this is the current energy, exit ungracefully
            if (!stretched) {

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
            else return numeric_limits<double>::infinity();
        } else energy += tempEnergy;
#ifdef SERIAL_CUDACROSSCHECK
        cudaDeviceSynchronize();
        resetdoublevariableCUDA<<<1,1,0, streamF>>>(gpu_Uvec);
        addvectorred3<<<bntaddvec2.at(2),bntaddvec2.at(3), bntaddvec2.at(3) * sizeof(double)
                                                                   , streamF>>>
        (CUDAcommon::cudavars.gpu_energyvec, gpu_params,
                gpu_Uvec);
        double cuda_energyvec[1];
        CUDAcommon::handleerror(cudaMemcpy(cuda_energyvec, gpu_Uvec, sizeof(double),
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
    chrono::duration<double> elapsed_run2(tend - tbegin);
    CUDAcommon::cudatime.TveccomputeE.push_back(elapsed_run2.count());
    CUDAcommon::cudatime.TcomputeE += elapsed_run2.count();
    CUDAcommon::cudatime.TcomputeEiter += elapsed_run2.count();
    tbegin = chrono::high_resolution_clock::now();
#endif
//    std::cout<<"CUDA energy total nint "<<CUDAcommon::cudavars.offset_E<<endl;
    /*vector<double> ones;
    for(int i = 0;i<8192;i++)
        ones.push_back(1.0);
    CUDAcommon::handleerror(cudaMemcpyAsync(CUDAcommon::cudavars.gpu_energyvec, ones
                                                        .data() ,
                                                bntaddvec2.at(0) * sizeof
                                                        (double),
                                                cudaMemcpyHostToDevice,streamF));*/
/*    CUDAcommon::handleerror(cudaMemsetAsync(CUDAcommon::cudavars.gpu_energyvec, 1,
                                            bntaddvec2.at(0) * sizeof
            (double),streamF));
    cudaDeviceSynchronize();*/
    resetdoublevariableCUDA<<<1,1,0, streamF>>>(gpu_Uvec);
    addvectorred3<<<bntaddvec2.at(2),bntaddvec2.at(3), bntaddvec2.at(3) * sizeof(double)
                    , streamF>>>(CUDAcommon::cudavars.gpu_energyvec, gpu_params,
                        gpu_Uvec);
    CUDAcommon::handleerror(cudaStreamSynchronize(streamF));
#ifdef DETAILEDOUTPUT_ENERGY
    cudaDeviceSynchronize();
    double cuda_energyvec[1];
    CUDAcommon::handleerror(cudaMemcpy(cuda_energyvec, gpu_Uvec, sizeof(double),
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
    chrono::duration<double> elapsed_run3(tend - tbegin);
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
    updateMembraneSharedData<stretched>();
    return energy;
}
template double ForceFieldManager::computeEnergy< false >(double *, bool) const;
template double ForceFieldManager::computeEnergy< true >(double *, bool) const;

void ForceFieldManager::computeForces(double *coord, double *f) {
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
    for (int i = 0; i < Bead::getDbData().forces.size_raw(); i++)
        f[i] = 0.0;
    //@}
#ifdef CUDATIMETRACK
    tend= chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed_run(tend - tbegin);
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
    chrono::duration<double> elapsed_run2(tend - tbegin);
    CUDAcommon::cudatime.TveccomputeF.push_back(elapsed_run2.count());
    CUDAcommon::cudatime.TcomputeF += elapsed_run2.count();
    tbegin = chrono::high_resolution_clock::now();
#endif

    // compute geometry
    updateGeometryValueWithDerivative();

    //recompute
//    double *F_i = new double[CGMethod::N];
    for (auto &ff : _forceFields) {
        ff->computeForces(coord, f);
#ifdef ALLSYNC
        cudaDeviceSynchronize();
#endif

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

    updateForceSharedData();
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
    copyForcesCUDA << < blocksnthreads[0], blocksnthreads[1], 0, stream >> >
                                                                 (f, fprev, gpu_nint);
    CUDAcommon::handleerror(cudaGetLastError(), "copyForcesCUDA", "ForceFieldManager.cu");
}

void ForceFieldManager::assignallforcemags() {

    for (auto &ff : _forceFields)
        ff->assignforcemags();
}

#endif
