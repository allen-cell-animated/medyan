
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

#include "CylinderExclVolume.h"

#include "CylinderExclVolRepulsion.h"

#include "Cylinder.h"
#include "Bead.h"

#include "MathFunctions.h"
#include "cross_check.h"
#include "nvToolsExt.h"

using namespace mathfunc;

template <class CVolumeInteractionType>
int CylinderExclVolume<CVolumeInteractionType>::numInteractions;

template <class CVolumeInteractionType>
void CylinderExclVolume<CVolumeInteractionType>::vectorize() {
    //count interactions
    nint = 0;

    for(auto ci : Cylinder::getCylinders()) {
        if(!ci->isFullLength()) continue;
        //do not calculate exvol for a non full length cylinder
        for(auto &cn : _neighborList->getNeighbors(ci))
        {
            if(!cn->isFullLength()||
               cn->getBranchingCylinder() == ci) continue;
            nint++;
        }
    }

    numInteractions = nint;
//    std::cout<<"NINT1 "<<nint<<endl;
    beadSet = new int[n * nint];
    krep = new double[nint];


    int nc = Cylinder::getCylinders().size();
    int i = 0;
    int Cumnc=0;
    for (i = 0; i < nc; i++) {
        auto ci = Cylinder::getCylinders()[i];
        if(!ci->isFullLength()) continue;

        int nn = _neighborList->getNeighbors(ci).size();
//        std::cout<<"Cylinder "<<i<<" "<<nn<<endl;
        for (int ni = 0; ni < nn; ni++) {

            auto cin = _neighborList->getNeighbors(ci)[ni];
            if(!cin->isFullLength()||
               cin->getBranchingCylinder() == ci) continue;
            beadSet[n * (Cumnc)] = ci->getFirstBead()->_dbIndex;
            beadSet[n * (Cumnc) + 1] = ci->getSecondBead()->_dbIndex;
            beadSet[n * (Cumnc) + 2] = cin->getFirstBead()->_dbIndex;
            beadSet[n * (Cumnc) + 3] = cin->getSecondBead()->_dbIndex;
            krep[Cumnc] = ci->getMCylinder()->getExVolConst();
            Cumnc++;
        }
    }
    //CUDA
#ifdef CUDAACCL
    nvtxRangePushA("CVFF");

//    cudaEvent_t start, stop;
//    CUDAcommon::handleerror(cudaEventCreate( &start));
//    CUDAcommon::handleerror(cudaEventCreate( &stop));
//    CUDAcommon::handleerror(cudaEventRecord( start, 0));

//    blocksnthreads.push_back(int(numInteractions/THREADSPERBLOCK + 1));
//    if(blocksnthreads[0]==1) blocksnthreads.push_back( numInteractions);
////    if(blocksnthreads[0]==1) blocksnthreads.push_back( 32*(int(numInteractions/32 +1)) );
//    else blocksnthreads.push_back(THREADSPERBLOCK);
    _FFType.optimalblocksnthreads(numInteractions);

    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_beadSet, n * numInteractions * sizeof(int)),"cuda data transfer",
                            "CylinderExclVolume.cu");
    CUDAcommon::handleerror(cudaMemcpy(gpu_beadSet, beadSet, n * numInteractions * sizeof(int), cudaMemcpyHostToDevice),
                            "cuda data transfer", "CylinderExclVolume.cu");
    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_krep, numInteractions * sizeof(double)),"cuda data transfer",
                            "CylinderExclVolume.cu");
    CUDAcommon::handleerror(cudaMemcpy(gpu_krep, krep, numInteractions * sizeof(double), cudaMemcpyHostToDevice),
                            "cuda data transfer", "CylinderExclVolume.cu");
    vector<int> params;
    params.push_back(int(n));
    params.push_back(numInteractions);
    //TODO make sure not using cudafree here is fine.
    if(gpu_params != NULL )
        CUDAcommon::handleerror(cudaFree(gpu_params),"cudaFree", "CylinderExclVolume.cu");
    if(nint > 0) {
        CUDAcommon::handleerror(cudaMalloc((void **) &gpu_params, 2 * sizeof(int)), "cuda data transfer",
                                "CylinderExclVolume.cu");
        CUDAcommon::handleerror(cudaMemcpy(gpu_params, params.data(), 2 * sizeof(int), cudaMemcpyHostToDevice),
                                "cuda data transfer", "CylinderExclVolume.cu");
    }
    nvtxRangePop();
#endif
    //
}


template <class CVolumeInteractionType>
void CylinderExclVolume<CVolumeInteractionType>::deallocate() {

    delete beadSet;
    delete krep;
#ifdef CUDAACCL
    if(nint > 0) {
        _FFType.deallocate();
        CUDAcommon::handleerror(cudaFree(gpu_beadSet), "cudaFree", "CylinderExclVolume.cu");
        CUDAcommon::handleerror(cudaFree(gpu_krep), "cudaFree", "CylinderExclVolume.cu");
        CUDAcommon::handleerror(cudaFree(gpu_params), "cudaFree", "CylinderExclVolume.cu");
        gpu_beadSet = NULL;
        gpu_krep = NULL;
        gpu_params = NULL;
    }
#endif
}


template <class CVolumeInteractionType>
double CylinderExclVolume<CVolumeInteractionType>::computeEnergy(double *coord, double *f, double d) {

    double U_i[1];
    double U_ii=NULL;
    double *gU_i;
#ifdef CUDAACCL
    //has to be changed to accomodate aux force
    double * gpu_coord=CUDAcommon::getCUDAvars().gpu_coord;
    double * gpu_force = CUDAcommon::getCUDAvars().gpu_force;
    double * gpu_d = CUDAcommon::getCUDAvars().gpu_lambda;
    nvtxRangePushA("CCEE");
//    if(d == 0.0){
//        gU_i=_FFType.energy(gpu_coord, gpu_force, gpu_beadSet, gpu_krep, gpu_params);
//    }
//    else{
        gU_i=_FFType.energy(gpu_coord, gpu_force, gpu_beadSet, gpu_krep, gpu_d, gpu_params);
//    }
    nvtxRangePop();
#else
    nvtxRangePushA("SCEE");
    if (d == 0.0)
        U_ii = _FFType.energy(coord, f, beadSet, krep);
    else
        U_ii = _FFType.energy(coord, f, beadSet, krep, d);
//        U_i = _FFType.energy(coord, f, beadSet, krep, d);
    nvtxRangePop();
#endif
    return U_ii;
}

template <class CVolumeInteractionType>
void CylinderExclVolume<CVolumeInteractionType>::computeForces(double *coord, double *f) {

#ifdef CUDAACCL
    //has to be changed to accomodate aux force
    double * gpu_coord=CUDAcommon::getCUDAvars().gpu_coord;
    //
    double * gpu_force;
//    std::cout<<"ex  vol cda"<<endl;
    if(cross_checkclass::Aux) {
        nvtxRangePushA("CCFE");
        gpu_force = CUDAcommon::getCUDAvars().gpu_forceAux;
        _FFType.forces(gpu_coord, gpu_force, gpu_beadSet, gpu_krep, gpu_params);
        nvtxRangePop();
    }
    else {
        nvtxRangePushA("CCFE");
        gpu_force = CUDAcommon::getCUDAvars().gpu_force;
        _FFType.forces(gpu_coord, gpu_force, gpu_beadSet, gpu_krep, gpu_params);
        nvtxRangePop();
    }
#else
    nvtxRangePushA("SCFE");
//    std::cout<<"x vol srl"<<endl;
    _FFType.forces(coord, f, beadSet, krep);
    nvtxRangePop();
#endif
}

///Template specializations
template double CylinderExclVolume<CylinderExclVolRepulsion>::computeEnergy(double *coord, double *f, double d);
template void CylinderExclVolume<CylinderExclVolRepulsion>::computeForces(double *coord, double *f);
template void CylinderExclVolume<CylinderExclVolRepulsion>::vectorize();
template void CylinderExclVolume<CylinderExclVolRepulsion>::deallocate();


