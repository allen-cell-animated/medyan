
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

#include "BranchingBending.h"

#include "BranchingBendingCosine.h"

#include "BranchingPoint.h"
#include "Cylinder.h"
#include "Bead.h"
#include "nvToolsExt.h"
#include "cross_check.h"

template <class BBendingInteractionType>
void BranchingBending<BBendingInteractionType>::vectorize() {

    beadSet = new int[n * BranchingPoint::getBranchingPoints().size()];
    kbend = new double[BranchingPoint::getBranchingPoints().size()];
    eqt = new double[BranchingPoint::getBranchingPoints().size()];

    int i = 0;

    for (auto b: BranchingPoint::getBranchingPoints()) {

        beadSet[n * i] = b->getFirstCylinder()->getFirstBead()->_dbIndex;
        beadSet[n * i + 1] = b->getFirstCylinder()->getSecondBead()->_dbIndex;
        beadSet[n * i + 2] = b->getSecondCylinder()->getFirstBead()->_dbIndex;
        beadSet[n * i + 3] = b->getSecondCylinder()->getSecondBead()->_dbIndex;

        kbend[i] = b->getMBranchingPoint()->getStretchingConstant();
        eqt[i] = b->getMBranchingPoint()->getEqTheta();

        i++;
    }
    //CUDA
#ifdef CUDAACCL
    F_i = new double [3 * Bead::getBeads().size()];
    nvtxRangePushA("CVFF");
    int numInteractions = BranchingPoint::getBranchingPoints().size();
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
    nvtxRangePop();
#endif
}

template<class BBendingInteractionType>
void BranchingBending<BBendingInteractionType>::deallocate() {

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



template <class BBendingInteractionType>
double BranchingBending<BBendingInteractionType>::computeEnergy(double *coord, double *f, double d) {

    double U_i[1], U_ii;
    double* gU_i;
    U_ii = NULL;
#ifdef CUDAACCL
    //has to be changed to accomodate aux force
    double * gpu_coord=CUDAcommon::getCUDAvars().gpu_coord;
    double * gpu_force=CUDAcommon::getCUDAvars().gpu_force;
    double * gpu_d = CUDAcommon::getCUDAvars().gpu_lambda;
    nvtxRangePushA("CCEBB");

    if(d == 0.0){
        nvtxRangePushA("CCEBB1");
        gU_i=_FFType.energy(gpu_coord, gpu_force, gpu_beadSet, gpu_kbend, gpu_eqt, gpu_params);
        nvtxRangePop();
    }
    else{
        nvtxRangePushA("CCEBB2");
        gU_i=_FFType.energy(gpu_coord, gpu_force, gpu_beadSet, gpu_kbend, gpu_eqt, gpu_d,
                            gpu_params);
        nvtxRangePop();
    }
    nvtxRangePop();
#else
    nvtxRangePushA("SCEBB");

    if (d == 0.0)
        U_ii = _FFType.energy(coord, f, beadSet, kbend, eqt);
    else
        U_ii = _FFType.energy(coord, f, beadSet, kbend, eqt, d);
    nvtxRangePop();
#endif
    return U_ii;
}

template <class BBendingInteractionType>
void BranchingBending<BBendingInteractionType>::computeForces(double *coord, double *f) {
#ifdef CUDAACCL
    //has to be changed to accomodate aux force
    double * gpu_coord=CUDAcommon::getCUDAvars().gpu_coord;

    double * gpu_force;

    if(cross_checkclass::Aux){
        nvtxRangePushA("CCFBB");

        gpu_force=CUDAcommon::getCUDAvars().gpu_forceAux;
        _FFType.forces(gpu_coord, gpu_force, gpu_beadSet, gpu_kbend, gpu_eqt, gpu_params);
        nvtxRangePop();
    }
    else {
        nvtxRangePushA("CCFBB");

        gpu_force = CUDAcommon::getCUDAvars().gpu_force;
        _FFType.forces(gpu_coord, gpu_force, gpu_beadSet, gpu_kbend, gpu_eqt, gpu_params);
        nvtxRangePop();
    }
#else
    nvtxRangePushA("SCFBB");

    _FFType.forces(coord, f, beadSet, kbend, eqt);
    nvtxRangePop();
#endif
}

///Template specializations
template double BranchingBending<BranchingBendingCosine>::computeEnergy(double *coord, double *f, double d);
template void BranchingBending<BranchingBendingCosine>::computeForces(double *coord, double *f);
template void BranchingBending<BranchingBendingCosine>::vectorize();
template void BranchingBending<BranchingBendingCosine>::deallocate();
