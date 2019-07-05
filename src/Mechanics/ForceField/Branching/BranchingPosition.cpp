
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v4.0
//
//  Copyright (2015-2018)  Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the MEDYAN web page for more information:
//  http://www.medyan.org
//------------------------------------------------------------------

#include "BranchingPosition.h"

#include "BranchingPositionCosine.h"

#include "BranchingPoint.h"
#include "Cylinder.h"
#include "Bead.h"

#ifdef CUDAACCL
#include "nvToolsExt.h"
#include "cross_check.h"
#endif
#include "Mechanics/CUDAcommon.h"

template <class BPositionInteractionType>
void BranchingPosition<BPositionInteractionType>::vectorize() {
    CUDAcommon::tmin.numinteractions[7] += BranchingPoint::getBranchingPoints().size();
    beadSet = new int[n * BranchingPoint::getBranchingPoints().size()];
    kpos = new floatingpoint[BranchingPoint::getBranchingPoints().size()];
    pos = new floatingpoint[BranchingPoint::getBranchingPoints().size()];

    int i = 0;

    for (auto b: BranchingPoint::getBranchingPoints()) {

        beadSet[n * i] = b->getFirstCylinder()->getFirstBead()->getStableIndex();
        beadSet[n * i + 1] = b->getFirstCylinder()->getSecondBead()->getStableIndex();
        beadSet[n * i + 2] = b->getSecondCylinder()->getFirstBead()->getStableIndex();

        kpos[i] = b->getMBranchingPoint()->getPositionConstant();
        pos[i] = b->getPosition();

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

    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_kpos, numInteractions * sizeof(floatingpoint)));
    CUDAcommon::handleerror(cudaMemcpy(gpu_kpos, kpos, numInteractions * sizeof(floatingpoint), cudaMemcpyHostToDevice));

    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_pos, numInteractions * sizeof(floatingpoint)));
    CUDAcommon::handleerror(cudaMemcpy(gpu_pos, pos, numInteractions * sizeof(floatingpoint), cudaMemcpyHostToDevice));

    vector<int> params;
    params.push_back(int(n));
    params.push_back(numInteractions);
    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_params, 2 * sizeof(int)));
    CUDAcommon::handleerror(cudaMemcpy(gpu_params, params.data(), 2 * sizeof(int), cudaMemcpyHostToDevice));
#endif

    //
}

template<class BPositionInteractionType>
void BranchingPosition<BPositionInteractionType>::deallocate() {
    delete [] beadSet;
    delete [] kpos;
    delete [] pos;
#ifdef CUDAACCL
    _FFType.deallocate();
    CUDAcommon::handleerror(cudaFree(gpu_beadSet));
    CUDAcommon::handleerror(cudaFree(gpu_kpos));
    CUDAcommon::handleerror(cudaFree(gpu_pos));
    CUDAcommon::handleerror(cudaFree(gpu_params));
#endif
}

template <class BPositionInteractionType>
floatingpoint BranchingPosition<BPositionInteractionType>::computeEnergy(floatingpoint *coord) {

    floatingpoint U_ii=(floatingpoint) 0.0;

#ifdef CUDAACCL
    //has to be changed to accomodate aux force
    floatingpoint * gpu_coord=CUDAcommon::getCUDAvars().gpu_coord;
    floatingpoint * gpu_force=CUDAcommon::getCUDAvars().gpu_force;
    floatingpoint * gpu_d = CUDAcommon::getCUDAvars().gpu_lambda;

//    if(d == 0.0){
//        gU_i=_FFType.energy(gpu_coord, gpu_force, gpu_beadSet, gpu_kpos, gpu_pos, gpu_params);
//
//    }
//    else{
        gU_i=_FFType.energy(gpu_coord, gpu_force, gpu_beadSet, gpu_kpos, gpu_pos, gpu_d,
                            gpu_params);
//    }

#endif
#ifdef SERIAL

    U_ii = _FFType.energy(coord, beadSet, kpos, pos);

#endif
#if defined(SERIAL_CUDACROSSCHECK) && defined(DETAILEDOUTPUT_ENERGY)
	floatingpoint U_i[1];
	floatingpoint* gU_i;
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

template <class BPositionInteractionType>
void BranchingPosition<BPositionInteractionType>::computeForces(floatingpoint *coord, floatingpoint *f) {
#ifdef CUDAACCL
    //has to be changed to accomodate aux force
    floatingpoint * gpu_coord=CUDAcommon::getCUDAvars().gpu_coord;

    floatingpoint * gpu_force;


    if(cross_checkclass::Aux){


        gpu_force=CUDAcommon::getCUDAvars().gpu_forceAux;
        _FFType.forces(gpu_coord, gpu_force, gpu_beadSet, gpu_kpos, gpu_pos, gpu_params);

    }
    else {

        gpu_force = CUDAcommon::getCUDAvars().gpu_force;
        _FFType.forces(gpu_coord, gpu_force, gpu_beadSet, gpu_kpos, gpu_pos, gpu_params);
    }
#endif
#ifdef SERIAL

    _FFType.forces(coord, f, beadSet, kpos, pos);

#endif
}


///Template specializations
template floatingpoint BranchingPosition<BranchingPositionCosine>::computeEnergy(floatingpoint *coord);
template void BranchingPosition<BranchingPositionCosine>::computeForces(floatingpoint *coord, floatingpoint *f);
template void BranchingPosition<BranchingPositionCosine>::vectorize();
template void BranchingPosition<BranchingPositionCosine>::deallocate();
