
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
#include "FilamentBending.h"

#include "FilamentBendingHarmonic.h"
#include "FilamentBendingCosine.h"

#include "Filament.h"
#include "Cylinder.h"
#include "Bead.h"
#include "nvToolsExt.h"
#include "cross_check.h"
template <class FBendingInteractionType>
void FilamentBending<FBendingInteractionType>::vectorize() {

    int numInteractions = Bead::getBeads().size() - 2 * Filament::getFilaments().size();

    beadSet = new int[n * numInteractions];
    kbend = new double[numInteractions];
    eqt = new double[numInteractions];

    int i = 0;

    for (auto f: Filament::getFilaments()) {

        if (f->getCylinderVector().size() > 1){

            for (auto it = f->getCylinderVector().begin()+1;
                 it != f->getCylinderVector().end(); it++){

                auto it2 = it - 1;
                beadSet[n * i] = (*it2)->getFirstBead()->_dbIndex;
                beadSet[n * i + 1] = (*it)->getFirstBead()->_dbIndex;;
                beadSet[n * i + 2] = (*it)->getSecondBead()->_dbIndex;;

                kbend[i] = (*it)->getMCylinder()->getBendingConst();
                eqt[i]  = (*it)->getMCylinder()->getEqTheta();

                i++;
            }
        }
    }

    //CUDA
#ifdef CUDAACCL
//    F_i = new double[3 * Bead::getBeads().size()];
//    nvtxRangePushA("CVFF");
    _FFType.optimalblocksnthreads(numInteractions);

    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_beadSet, n * numInteractions * sizeof(int)));
    CUDAcommon::handleerror(cudaMemcpy(gpu_beadSet, beadSet, n * numInteractions * sizeof(int),
                                       cudaMemcpyHostToDevice));

    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_kbend, numInteractions * sizeof(double)));
    CUDAcommon::handleerror(cudaMemcpy(gpu_kbend, kbend, numInteractions * sizeof(double), cudaMemcpyHostToDevice));

    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_eqt, numInteractions * sizeof(double)));
    CUDAcommon::handleerror(cudaMemcpy(gpu_eqt, eqt, numInteractions * sizeof(double), cudaMemcpyHostToDevice));

    vector<int> params;
    params.push_back(int(n));
    params.push_back(numInteractions);
    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_params, 2 * sizeof(int)));
    CUDAcommon::handleerror(cudaMemcpy(gpu_params, params.data(), 2 * sizeof(int), cudaMemcpyHostToDevice));
//    nvtxRangePop();
#endif
}

template<class FBendingInteractionType>
void FilamentBending<FBendingInteractionType>::deallocate() {

    delete beadSet;
    delete kbend;
    delete eqt;
#ifdef CUDAACCL
    _FFType.deallocate();
    CUDAcommon::handleerror(cudaFree(gpu_beadSet));
    CUDAcommon::handleerror(cudaFree(gpu_kbend));
    CUDAcommon::handleerror(cudaFree(gpu_eqt));
    CUDAcommon::handleerror(cudaFree(gpu_params));
#endif
}


template <class FBendingInteractionType>
double FilamentBending<FBendingInteractionType>::computeEnergy(double *coord, double *f, double d){

    double U_i[1], U_ii;
    double* gU_i;
    U_ii = NULL;
#ifdef CUDAACCL
    //has to be changed to accomodate aux force
    double * gpu_coord=CUDAcommon::getCUDAvars().gpu_coord;
    double * gpu_force=CUDAcommon::getCUDAvars().gpu_force;
    double * gpu_d = CUDAcommon::getCUDAvars().gpu_lambda;
//    nvtxRangePushA("CCEFB");

//    if(d == 0.0){
//        gU_i=_FFType.energy(gpu_coord, gpu_force, gpu_beadSet, gpu_kbend, gpu_eqt, gpu_params);
//
//    }
//    else{
        gU_i=_FFType.energy(gpu_coord, gpu_force, gpu_beadSet, gpu_kbend, gpu_eqt, gpu_d,
                            gpu_params);
//    }
//    nvtxRangePop();
#endif
#ifdef SERIAL
//    nvtxRangePushA("SCEFB");
    if (d == 0.0)
        U_ii = _FFType.energy(coord, f, beadSet, kbend, eqt);
    else
        U_ii= _FFType.energy(coord, f, beadSet, kbend, eqt, d);
//    nvtxRangePop();
#endif
    return U_ii;

}

template <class FBendingInteractionType>
void FilamentBending<FBendingInteractionType>::computeForces(double *coord, double *f) {
#ifdef CUDAACCL
    //has to be changed to accomodate aux force
    double * gpu_coord=CUDAcommon::getCUDAvars().gpu_coord;

    double * gpu_force;


    if(cross_checkclass::Aux){
//        nvtxRangePushA("CCFFB");

        gpu_force=CUDAcommon::getCUDAvars().gpu_forceAux;
        _FFType.forces(gpu_coord, gpu_force, gpu_beadSet, gpu_kbend, gpu_eqt, gpu_params);
//        nvtxRangePop();
    }
    else {
//        nvtxRangePushA("CCFFB");

        gpu_force = CUDAcommon::getCUDAvars().gpu_force;
        _FFType.forces(gpu_coord, gpu_force, gpu_beadSet, gpu_kbend, gpu_eqt, gpu_params);
//        nvtxRangePop();
    }
#endif
#ifdef SERIAL
//    nvtxRangePushA("SCFFB");

    _FFType.forces(coord, f, beadSet, kbend, eqt);
//    nvtxRangePop();
#endif

}

///Template specializations
template double FilamentBending<FilamentBendingHarmonic>::computeEnergy(double *coord, double *f, double d);
template void FilamentBending<FilamentBendingHarmonic>::computeForces(double *coord, double *f);
template void FilamentBending<FilamentBendingHarmonic>::vectorize();
template void FilamentBending<FilamentBendingHarmonic>::deallocate();


template double FilamentBending<FilamentBendingCosine>::computeEnergy(double *coord, double *f, double d);
template void FilamentBending<FilamentBendingCosine>::computeForces(double *coord, double *f);
template void FilamentBending<FilamentBendingCosine>::vectorize();
template void FilamentBending<FilamentBendingCosine>::deallocate();
