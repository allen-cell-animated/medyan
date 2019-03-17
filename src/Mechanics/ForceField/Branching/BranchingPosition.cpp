
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

#include "BranchingPosition.h"

#include "BranchingPositionCosine.h"

#include "BranchingPoint.h"
#include "Cylinder.h"
#include "Bead.h"

#ifdef CUDAACCL
#include "nvToolsExt.h"
#include "cross_check.h"
#endif

template <class BPositionInteractionType>
void BranchingPosition<BPositionInteractionType>::vectorize() {

    beadSet = new int[n * BranchingPoint::getBranchingPoints().size()];
    kpos = new double[BranchingPoint::getBranchingPoints().size()];
    pos = new double[BranchingPoint::getBranchingPoints().size()];

    int i = 0;

    for (auto b: BranchingPoint::getBranchingPoints()) {

        beadSet[n * i] = b->getFirstCylinder()->getFirstBead()->getDbIndex();
        beadSet[n * i + 1] = b->getFirstCylinder()->getSecondBead()->getDbIndex();
        beadSet[n * i + 2] = b->getSecondCylinder()->getFirstBead()->getDbIndex();

        kpos[i] = b->getMBranchingPoint()->getPositionConstant();
        pos[i] = b->getPosition();

        i++;
    }
    //CUDA
#ifdef CUDAACCL
//    F_i = new double [3 * Bead::getBeads().size()];
    int numInteractions = BranchingPoint::getBranchingPoints().size();
    _FFType.optimalblocksnthreads(numInteractions);

    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_beadSet, n * numInteractions * sizeof(int)));
    CUDAcommon::handleerror(cudaMemcpy(gpu_beadSet, beadSet, n * numInteractions * sizeof(int),
                                       cudaMemcpyHostToDevice));

    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_kpos, numInteractions * sizeof(double)));
    CUDAcommon::handleerror(cudaMemcpy(gpu_kpos, kpos, numInteractions * sizeof(double), cudaMemcpyHostToDevice));

    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_pos, numInteractions * sizeof(double)));
    CUDAcommon::handleerror(cudaMemcpy(gpu_pos, pos, numInteractions * sizeof(double), cudaMemcpyHostToDevice));

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
double BranchingPosition<BPositionInteractionType>::computeEnergy(double *coord, double *f, double d) {

    double U_i[1], U_ii;
    double* gU_i;
    U_ii = 0.0;
#ifdef CUDAACCL
    //has to be changed to accomodate aux force
    double * gpu_coord=CUDAcommon::getCUDAvars().gpu_coord;
    double * gpu_force=CUDAcommon::getCUDAvars().gpu_force;
    double * gpu_d = CUDAcommon::getCUDAvars().gpu_lambda;

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
    if (d == 0.0)
        U_ii = _FFType.energy(coord, f, beadSet, kpos, pos);
    else
        U_ii = _FFType.energy(coord, f, beadSet, kpos, pos, d);
#endif
#if defined(SERIAL_CUDACROSSCHECK) && defined(DETAILEDOUTPUT_ENERGY)
    CUDAcommon::handleerror(cudaDeviceSynchronize(),"ForceField", "ForceField");
    double cuda_energy[1];
    if(gU_i == NULL)
        cuda_energy[0] = 0.0;
    else {
        CUDAcommon::handleerror(cudaMemcpy(cuda_energy, gU_i, sizeof(double),
                                           cudaMemcpyDeviceToHost));
    }
    std::cout<<getName()<<" Serial Energy "<<U_ii<<" Cuda Energy "<<cuda_energy[0]<<endl;
#endif
    return U_ii;
}

template <class BPositionInteractionType>
void BranchingPosition<BPositionInteractionType>::computeForces(double *coord, double *f) {
#ifdef CUDAACCL
    //has to be changed to accomodate aux force
    double * gpu_coord=CUDAcommon::getCUDAvars().gpu_coord;

    double * gpu_force;


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
template double BranchingPosition<BranchingPositionCosine>::computeEnergy(double *coord, double *f, double d);
template void BranchingPosition<BranchingPositionCosine>::computeForces(double *coord, double *f);
template void BranchingPosition<BranchingPositionCosine>::vectorize();
template void BranchingPosition<BranchingPositionCosine>::deallocate();
