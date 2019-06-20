
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
#include "FilamentBending.h"

#include "FilamentBendingHarmonic.h"
#include "FilamentBendingCosine.h"

#include "Filament.h"
#include "Cylinder.h"
#include "Bead.h"
#include "Bubble.h"
#include "CGMethod.h"
#ifdef CUDAACCL
#include "nvToolsExt.h"
#endif
#include "cross_check.h"
template <class FBendingInteractionType>
void FilamentBending<FBendingInteractionType>::vectorize() {

    int numInteractions = Bead::getBeads().size() - 2 * Filament::getFilaments().size() - Bubble::numBubbles();

    beadSet = new int[n * numInteractions];
    kbend = new floatingpoint[numInteractions];
    eqt = new floatingpoint[numInteractions];

    int i = 0;

    for (auto f: Filament::getFilaments()) {

        if (f->getCylinderVector().size() > 1){

            for (auto it = f->getCylinderVector().begin()+1;
                 it != f->getCylinderVector().end(); it++){

                auto it2 = it - 1;
                beadSet[n * i] = (*it2)->getFirstBead()->getIndex();
                beadSet[n * i + 1] = (*it)->getFirstBead()->getIndex();
                beadSet[n * i + 2] = (*it)->getSecondBead()->getIndex();
//                std::cout<<f->getCylinderVector().size()<<" "<<(*it2)->getFirstBead()
//                        ->getIndex()<<" "<<(*it)->getFirstBead()
//                        ->getIndex()<<" "<<(*it)->getSecondBead()->getIndex()<<endl;
                kbend[i] = (*it)->getMCylinder()->getBendingConst();
                eqt[i]  = (*it)->getMCylinder()->getEqTheta();

                i++;
            }
        }
    }

    //CUDA
#ifdef CUDAACCL
#ifdef CUDATIMETRACK
    chrono::high_resolution_clock::time_point tbegin, tend;
    tbegin = chrono::high_resolution_clock::now();
#endif
    //CUDA stream create
    if(stream == NULL || !(CUDAcommon::getCUDAvars().conservestreams))
        CUDAcommon::handleerror(cudaStreamCreate(&stream));

//    F_i = new floatingpoint[3 * Bead::getBeads().size()];

    _FFType.optimalblocksnthreads(numInteractions, stream);

    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_beadSet, n * numInteractions * sizeof(int)));
    CUDAcommon::handleerror(cudaMemcpyAsync(gpu_beadSet, beadSet, n * numInteractions *
                                                sizeof(int),
                                       cudaMemcpyHostToDevice, stream));

    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_kbend, numInteractions * sizeof(floatingpoint)));
    CUDAcommon::handleerror(cudaMemcpyAsync(gpu_kbend, kbend, numInteractions * sizeof
                            (floatingpoint), cudaMemcpyHostToDevice, stream));

    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_eqt, numInteractions * sizeof(floatingpoint)));
    CUDAcommon::handleerror(cudaMemcpyAsync(gpu_eqt, eqt, numInteractions * sizeof(floatingpoint),
                                        cudaMemcpyHostToDevice, stream));

    vector<int> params;
    params.push_back(int(n));
    params.push_back(numInteractions);
    params.push_back(CUDAcommon::cudavars.offset_E);
    //set offset
    CUDAcommon::cudavars.offset_E += numInteractions;
//    std::cout<<"offset "<<getName()<<" "<<CUDAcommon::cudavars.offset_E<<endl;
//    std::cout<<"offset "<<getName()<<" "<<CUDAcommon::cudavars.offset_E<<endl;
    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_params, 3 * sizeof(int)));
    CUDAcommon::handleerror(cudaMemcpyAsync(gpu_params, params.data(), 3 * sizeof(int),
                                       cudaMemcpyHostToDevice, stream));

#ifdef CUDATIMETRACK
//    CUDAcommon::handleerror(cudaDeviceSynchronize(),"FilamentBending.cu",
//                            "vectorizeFF.cu");
    tend= chrono::high_resolution_clock::now();
    chrono::duration<floatingpoint> elapsed_run(tend - tbegin);
    CUDAcommon::cudatime.TvecvectorizeFF.push_back(elapsed_run.count());
    CUDAcommon::cudatime.TvectorizeFF += elapsed_run.count();
#endif
#endif
}

template<class FBendingInteractionType>
void FilamentBending<FBendingInteractionType>::deallocate() {

    delete [] beadSet;
    delete [] kbend;
    delete [] eqt;
#ifdef CUDAACCL
    if(!(CUDAcommon::getCUDAvars().conservestreams))
        CUDAcommon::handleerror(cudaStreamDestroy(stream));
    _FFType.deallocate();
    CUDAcommon::handleerror(cudaFree(gpu_beadSet));
    CUDAcommon::handleerror(cudaFree(gpu_kbend));
    CUDAcommon::handleerror(cudaFree(gpu_eqt));
    CUDAcommon::handleerror(cudaFree(gpu_params));
#endif
}


template <class FBendingInteractionType>
floatingpoint FilamentBending<FBendingInteractionType>::computeEnergy(floatingpoint *coord, floatingpoint *f, floatingpoint d){

    floatingpoint U_i[1], U_ii=0.0;
    floatingpoint* gU_i;
    U_ii = 0.0;
#ifdef CUDATIMETRACK
    chrono::high_resolution_clock::time_point tbegin, tend;
#endif
#ifdef CUDAACCL
#ifdef CUDATIMETRACK

    tbegin = chrono::high_resolution_clock::now();
#endif

    //has to be changed to accomodate aux force
    floatingpoint * gpu_coord=CUDAcommon::getCUDAvars().gpu_coord;
    floatingpoint * gpu_force=CUDAcommon::getCUDAvars().gpu_force;
    floatingpoint * gpu_d = CUDAcommon::getCUDAvars().gpu_lambda;


//    if(d == 0.0){
//        gU_i=_FFType.energy(gpu_coord, gpu_force, gpu_beadSet, gpu_kbend, gpu_eqt, gpu_params);
//
//    }
//    else{
        gU_i=_FFType.energy(gpu_coord, gpu_force, gpu_beadSet, gpu_kbend, gpu_eqt, gpu_d,
                            gpu_params);
//    }
#ifdef CUDATIMETRACK
//    CUDAcommon::handleerror(cudaDeviceSynchronize(),"FilamentBending.cu",
//                            "computeEnergy");
    tend= chrono::high_resolution_clock::now();
    chrono::duration<floatingpoint> elapsed_run(tend - tbegin);
    CUDAcommon::cudatime.TveccomputeE.push_back(elapsed_run.count());
    CUDAcommon::cudatime.TcomputeE += elapsed_run.count();
    CUDAcommon::cudatime.TcomputeEiter += elapsed_run.count();
#endif
#endif
#ifdef SERIAL
#ifdef CUDATIMETRACK
    tbegin = chrono::high_resolution_clock::now();
#endif

    if (d == 0.0)
        U_ii = _FFType.energy(coord, f, beadSet, kbend, eqt);
    else
        U_ii= _FFType.energy(coord, f, beadSet, kbend, eqt, d);

#ifdef CUDATIMETRACK
    tend= chrono::high_resolution_clock::now();
    chrono::duration<floatingpoint> elapsed_runs(tend - tbegin);
    CUDAcommon::serltime.TveccomputeE.push_back(elapsed_runs.count());
    CUDAcommon::serltime.TcomputeE += elapsed_runs.count();
    CUDAcommon::serltime.TcomputeEiter += elapsed_runs.count();
#endif
#endif
    return U_ii;

}

template <class FBendingInteractionType>
void FilamentBending<FBendingInteractionType>::computeForces(floatingpoint *coord, floatingpoint *f) {
#ifdef CUDATIMETRACK
    chrono::high_resolution_clock::time_point tbegin, tend;
    tbegin = chrono::high_resolution_clock::now();
#endif
#ifdef CUDAACCL
    //has to be changed to accomodate aux force
    floatingpoint * gpu_coord=CUDAcommon::getCUDAvars().gpu_coord;
    floatingpoint * gpu_force;
    if(cross_checkclass::Aux){
        gpu_force=CUDAcommon::getCUDAvars().gpu_forceAux;
        _FFType.forces(gpu_coord, gpu_force, gpu_beadSet, gpu_kbend, gpu_eqt, gpu_params);
    }
    else {
        gpu_force = CUDAcommon::getCUDAvars().gpu_force;
        _FFType.forces(gpu_coord, gpu_force, gpu_beadSet, gpu_kbend, gpu_eqt, gpu_params);
    }
#endif
#ifdef CUDATIMETRACK
    tend= chrono::high_resolution_clock::now();
    chrono::duration<floatingpoint> elapsed_run(tend - tbegin);
    CUDAcommon::cudatime.TveccomputeF.push_back(elapsed_run.count());
    CUDAcommon::cudatime.TcomputeF += elapsed_run.count();
    tbegin = chrono::high_resolution_clock::now();
#endif
#ifdef SERIAL
    _FFType.forces(coord, f, beadSet, kbend, eqt);
#endif
#ifdef CUDATIMETRACK
    tend= chrono::high_resolution_clock::now();
    chrono::duration<floatingpoint> elapsed_runs(tend - tbegin);
    CUDAcommon::serltime.TveccomputeF.push_back(elapsed_runs.count());
    CUDAcommon::serltime.TcomputeF += elapsed_runs.count();
#endif
#ifdef DETAILEDOUTPUT
    floatingpoint maxF = 0.0;
    floatingpoint mag = 0.0;
    for(int i = 0; i < CGMethod::N/3; i++) {
        mag = 0.0;
        for(int j = 0; j < 3; j++)
            mag += f[3 * i + j]*f[3 * i + j];
        mag = sqrt(mag);
//        std::cout<<"SL "<<i<<" "<<mag*mag<<" "<<forceAux[3 * i]<<" "<<forceAux[3 * i + 1]<<" "<<forceAux[3 * i +
//                2]<<endl;
        if(mag > maxF) maxF = mag;
    }
    std::cout<<"max "<<getName()<<" "<<maxF<<endl;
#endif

}

///Template specializations
template floatingpoint FilamentBending<FilamentBendingHarmonic>::computeEnergy(floatingpoint *coord, floatingpoint *f, floatingpoint d);
template void FilamentBending<FilamentBendingHarmonic>::computeForces(floatingpoint *coord, floatingpoint *f);
template void FilamentBending<FilamentBendingHarmonic>::vectorize();
template void FilamentBending<FilamentBendingHarmonic>::deallocate();


template floatingpoint FilamentBending<FilamentBendingCosine>::computeEnergy(floatingpoint *coord, floatingpoint *f, floatingpoint d);
template void FilamentBending<FilamentBendingCosine>::computeForces(floatingpoint *coord, floatingpoint *f);
template void FilamentBending<FilamentBendingCosine>::vectorize();
template void FilamentBending<FilamentBendingCosine>::deallocate();
