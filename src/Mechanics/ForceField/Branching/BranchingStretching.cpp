
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

#include <CGMethod.h>
#include "BranchingStretching.h"

#include "BranchingStretchingHarmonic.h"

#include "BranchingPoint.h"
#include "Cylinder.h"
#include "Bead.h"
#include "cross_check.h"
#ifdef CUDAACCL
#include "nvToolsExt.h"
#endif

template <class BStretchingInteractionType>
void BranchingStretching<BStretchingInteractionType>::vectorize() {

    CUDAcommon::tmin.numinteractions[4] += BranchingPoint::getBranchingPoints().size();
    beadSet = new int[n * BranchingPoint::getBranchingPoints().size()];
    kstr = new floatingpoint[BranchingPoint::getBranchingPoints().size()];
    eql = new floatingpoint[BranchingPoint::getBranchingPoints().size()];
    pos = new floatingpoint[BranchingPoint::getBranchingPoints().size()];
    stretchforce = new floatingpoint[BranchingPoint::getBranchingPoints().size()];

    int i = 0;

    for (auto b: BranchingPoint::getBranchingPoints()) {

        beadSet[n * i] = b->getFirstCylinder()->getFirstBead()->getIndex() * 3 + si.bead;
        beadSet[n * i + 1] = b->getFirstCylinder()->getSecondBead()->getIndex() * 3 + si.bead;
        beadSet[n * i + 2] = b->getSecondCylinder()->getFirstBead()->getIndex() * 3 + si.bead;

        kstr[i] = b->getMBranchingPoint()->getStretchingConstant();
        eql[i] = b->getMBranchingPoint()->getEqLength();
        pos[i] = b->getPosition();
        stretchforce[i] = 0.0;
        i++;
    }
    //CUDA
#ifdef CUDAACCL
//    F_i = new floatingpoint[CGMethod::N];

    int numInteractions = BranchingPoint::getBranchingPoints().size();
    _FFType.optimalblocksnthreads(numInteractions);

    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_beadSet, n * numInteractions * sizeof(int)));
    CUDAcommon::handleerror(cudaMemcpy(gpu_beadSet, beadSet, n * numInteractions * sizeof(int),
                                       cudaMemcpyHostToDevice));

    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_kstr, numInteractions * sizeof(floatingpoint)));
    CUDAcommon::handleerror(cudaMemcpy(gpu_kstr, kstr, numInteractions * sizeof(floatingpoint), cudaMemcpyHostToDevice));

    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_eql, numInteractions * sizeof(floatingpoint)));
    CUDAcommon::handleerror(cudaMemcpy(gpu_eql, eql, numInteractions * sizeof(floatingpoint), cudaMemcpyHostToDevice));

    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_pos, numInteractions * sizeof(floatingpoint)));
    CUDAcommon::handleerror(cudaMemcpy(gpu_pos, pos, numInteractions * sizeof(floatingpoint), cudaMemcpyHostToDevice));

    vector<int> params;
    params.push_back(int(n));
    params.push_back(numInteractions);
    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_params, 2 * sizeof(int)));
    CUDAcommon::handleerror(cudaMemcpy(gpu_params, params.data(), 2 * sizeof(int), cudaMemcpyHostToDevice));
#endif
}

template<class BStretchingInteractionType>
void BranchingStretching<BStretchingInteractionType>::deallocate() {
    int i = 0;
    for(auto b:BranchingPoint::getBranchingPoints()){
        //Using += to ensure that the stretching forces are additive.
        b->getMBranchingPoint()->stretchForce += stretchforce[i];
        i++;
    }
    delete [] stretchforce;
    delete [] beadSet;
    delete [] kstr;
    delete [] eql;
    delete [] pos;
#ifdef CUDAACCL
    _FFType.deallocate();
    CUDAcommon::handleerror(cudaFree(gpu_beadSet));
    CUDAcommon::handleerror(cudaFree(gpu_kstr));
    CUDAcommon::handleerror(cudaFree(gpu_eql));
    CUDAcommon::handleerror(cudaFree(gpu_pos));
    CUDAcommon::handleerror(cudaFree(gpu_params));
#endif
}


template <class BStretchingInteractionType>
floatingpoint BranchingStretching<BStretchingInteractionType>::computeEnergy(floatingpoint *coord) {


    floatingpoint U_ii=(floatingpoint)0.0;

#ifdef CUDAACCL
    //has to be changed to accomodate aux force
    floatingpoint * gpu_coord=CUDAcommon::getCUDAvars().gpu_coord;
    floatingpoint * gpu_force=CUDAcommon::getCUDAvars().gpu_force;
    floatingpoint * gpu_d = CUDAcommon::getCUDAvars().gpu_lambda;

//    if(d == 0.0){
//        gU_i=_FFType.energy(gpu_coord, gpu_force, gpu_beadSet, gpu_kstr, gpu_eql, gpu_pos, gpu_params);
//
//    }
//    else{
        gU_i=_FFType.energy(gpu_coord, gpu_force, gpu_beadSet, gpu_kstr, gpu_eql, gpu_pos, gpu_d,
                            gpu_params);
//    }
#endif
#ifdef SERIAL

    U_ii = _FFType.energy(coord, beadSet, kstr, eql, pos);

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
        std::cout << getName()<<" "<<"Serial Energy " << U_ii << " Cuda Energy " <<
                                                                            cuda_energy[0] << endl;
#endif
    return U_ii;
}

template <class BStretchingInteractionType>
void BranchingStretching<BStretchingInteractionType>::computeForces(floatingpoint *coord, floatingpoint *f) {
#ifdef CUDAACCL
    //has to be changed to accomodate aux force
    floatingpoint * gpu_coord=CUDAcommon::getCUDAvars().gpu_coord;

    floatingpoint * gpu_force;
    if(cross_checkclass::Aux){

        gpu_force=CUDAcommon::getCUDAvars().gpu_forceAux;
        _FFType.forces(gpu_coord, gpu_force, gpu_beadSet, gpu_kstr, gpu_eql, gpu_pos, gpu_params);
    }
    else {

        gpu_force = CUDAcommon::getCUDAvars().gpu_force;
        _FFType.forces(gpu_coord, gpu_force, gpu_beadSet, gpu_kstr, gpu_eql, gpu_pos, gpu_params);
    }
#endif
#ifdef SERIAL

    _FFType.forces(coord, f, beadSet, kstr, eql, pos, stretchforce);

#endif
}
///Template specializations
template floatingpoint
BranchingStretching<BranchingStretchingHarmonic>::computeEnergy(floatingpoint *coord);
template void BranchingStretching<BranchingStretchingHarmonic>::computeForces(floatingpoint *coord, floatingpoint *f);
template void BranchingStretching<BranchingStretchingHarmonic>::vectorize();
template void BranchingStretching<BranchingStretchingHarmonic>::deallocate();