
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

#include "BranchingBending.h"

#include "BranchingBendingCosine.h"

#include "BranchingPoint.h"
#include "Cylinder.h"
#include "Bead.h"
#ifdef CUDAACCL
#include "nvToolsExt.h"
#endif
#include "cross_check.h"
#include "Mechanics/CUDAcommon.h"

template <class BBendingInteractionType>
void BranchingBending<BBendingInteractionType>::vectorize() {

    CUDAcommon::tmin.numinteractions[5] += BranchingPoint::getBranchingPoints().size();
    beadSet = new int[n * BranchingPoint::getBranchingPoints().size()];
    kbend = new floatingpoint[BranchingPoint::getBranchingPoints().size()];
    eqt = new floatingpoint[BranchingPoint::getBranchingPoints().size()];

    int i = 0;

    for (auto b: BranchingPoint::getBranchingPoints()) {

        beadSet[n * i] = b->getFirstCylinder()->getFirstBead()->getStableIndex();
        beadSet[n * i + 1] = b->getFirstCylinder()->getSecondBead()->getStableIndex();
        beadSet[n * i + 2] = b->getSecondCylinder()->getFirstBead()->getStableIndex();
        beadSet[n * i + 3] = b->getSecondCylinder()->getSecondBead()->getStableIndex();

        kbend[i] = b->getMBranchingPoint()->getStretchingConstant();
        eqt[i] = b->getMBranchingPoint()->getEqTheta();
        i++;
    }
    //CUDA
#ifdef CUDAACCL
//    F_i = new floatingpoint [3 * Bead::getBeads().size()];

    int numInteractions = BranchingPoint::getBranchingPoints().size();
    _FFType.optimalblocksnthreads(numInteractions);

    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_beadSet, n * numInteractions * sizeof(int)));
    CUDAcommon::handleerror(cudaMemcpy(gpu_beadSet, beadSet, n * numInteractions * sizeof(int),
                                       cudaMemcpyHostToDevice));

    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_kbend, numInteractions * sizeof(floatingpoint)));
    CUDAcommon::handleerror(cudaMemcpy(gpu_kbend, kbend, numInteractions * sizeof(floatingpoint), cudaMemcpyHostToDevice));

    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_eqt, numInteractions * sizeof(floatingpoint)));
    CUDAcommon::handleerror(cudaMemcpy(gpu_eqt, eqt, numInteractions * sizeof(floatingpoint), cudaMemcpyHostToDevice));

    vector<int> params;
    params.push_back(int(n));
    params.push_back(numInteractions);
    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_params, 2 * sizeof(int)));
    CUDAcommon::handleerror(cudaMemcpy(gpu_params, params.data(), 2 * sizeof(int), cudaMemcpyHostToDevice));

#endif
}

template<class BBendingInteractionType>
void BranchingBending<BBendingInteractionType>::deallocate() {
    delete [] beadSet;
    delete [] kbend;
    delete [] eqt;
#ifdef CUDAACCL
    _FFType.deallocate();
    CUDAcommon::handleerror(cudaFree(gpu_beadSet));
    CUDAcommon::handleerror(cudaFree(gpu_kbend));
    CUDAcommon::handleerror(cudaFree(gpu_eqt));
    CUDAcommon::handleerror(cudaFree(gpu_params));
#endif
}



template <class BBendingInteractionType>
floatingpoint BranchingBending<BBendingInteractionType>::computeEnergy(floatingpoint *coord) {

    floatingpoint U_ii=(floatingpoint)0.0;

#ifdef CUDAACCL
    floatingpoint* gU_i;
    //has to be changed to accomodate aux force
    floatingpoint * gpu_coord=CUDAcommon::getCUDAvars().gpu_coord;
    floatingpoint * gpu_force=CUDAcommon::getCUDAvars().gpu_force;
    floatingpoint * gpu_d = CUDAcommon::getCUDAvars().gpu_lambda;


//    if(d == 0.0){

//        gU_i=_FFType.energy(gpu_coord, gpu_force, gpu_beadSet, gpu_kbend, gpu_eqt, gpu_params);

//    }
//    else{

        gU_i=_FFType.energy(gpu_coord, gpu_force, gpu_beadSet, gpu_kbend, gpu_eqt, gpu_d,
                            gpu_params);

//    }

#endif
#ifdef SERIAL

    U_ii = _FFType.energy(coord, beadSet, kbend, eqt);

#endif
#if defined(SERIAL_CUDACROSSCHECK) && defined(DETAILEDOUTPUT_ENERGY)
    CUDAcommon::handleerror(cudaDeviceSynchronize(),"ForceField", "ForceField");
    floatingpoint cuda_energy[1];
    if(gU_i == NULL)
        cuda_energy[0] = 0.0;
    else {
        CUDAcommon::handleerror(cudaMemcpy(cuda_energy, gU_i, sizeof(floatingpoint),
                                           cudaMemcpyDeviceToHost));
    }
    std::cout<<getName()<<" Serial Energy "<<U_ii<<" Cuda Energy "<<cuda_energy[0]<<endl;
#endif
    return U_ii;
}

template <class BBendingInteractionType>
void BranchingBending<BBendingInteractionType>::computeForces(floatingpoint *coord, floatingpoint *f) {
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
#ifdef SERIAL


    _FFType.forces(coord, f, beadSet, kbend, eqt);

#endif
}

///Template specializations
template floatingpoint BranchingBending<BranchingBendingCosine>::computeEnergy(floatingpoint *coord);
template void BranchingBending<BranchingBendingCosine>::computeForces(floatingpoint *coord, floatingpoint *f);
template void BranchingBending<BranchingBendingCosine>::vectorize();
template void BranchingBending<BranchingBendingCosine>::deallocate();
