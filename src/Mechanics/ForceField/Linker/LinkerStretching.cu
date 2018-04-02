
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

#include "LinkerStretching.h"

#include "LinkerStretchingHarmonic.h"

#include "Cylinder.h"
#include "Linker.h"
#include "Bead.h"
#include "cross_check.h"
#include "nvToolsExt.h"
template <class LStretchingInteractionType>
void LinkerStretching<LStretchingInteractionType>::vectorize() {

    beadSet = new int[n * Linker::getLinkers().size()];
    kstr = new double[Linker::getLinkers().size()];
    eql = new double[Linker::getLinkers().size()];
    pos1 = new double[Linker::getLinkers().size()];
    pos2 = new double[Linker::getLinkers().size()];

    int i = 0;

    for (auto l: Linker::getLinkers()) {

        beadSet[n * i] = l->getFirstCylinder()->getFirstBead()->_dbIndex;
        beadSet[n * i + 1] = l->getFirstCylinder()->getSecondBead()->_dbIndex;
        beadSet[n * i + 2] = l->getSecondCylinder()->getFirstBead()->_dbIndex;
        beadSet[n * i + 3] = l->getSecondCylinder()->getSecondBead()->_dbIndex;

        kstr[i] = l->getMLinker()->getStretchingConstant();
        eql[i] = l->getMLinker()->getEqLength();
        pos1[i] = l->getFirstPosition();
        pos2[i] = l->getSecondPosition();

        i++;
    }

    //CUDA
#ifdef CUDAACCL
//    F_i = new double[3 * Bead::getBeads().size()];
//    cudaEvent_t start, stop;
//    CUDAcommon::handleerror(cudaEventCreate( &start));
//    CUDAcommon::handleerror(cudaEventCreate( &stop));
//    CUDAcommon::handleerror(cudaEventRecord( start, 0));
    nvtxRangePushA("CVFF");

    int numInteractions =Linker::getLinkers().size();
    _FFType.optimalblocksnthreads(numInteractions);
    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_beadSet, n * numInteractions * sizeof(int)),"cuda data transfer",
                                       "LinkerStretching.cu");
    CUDAcommon::handleerror(cudaMemcpy(gpu_beadSet, beadSet, n * numInteractions * sizeof(int),
                                       cudaMemcpyHostToDevice),"cuda data transfer", "LinkerStretching.cu");

    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_kstr, numInteractions * sizeof(double)),"cuda data transfer",
                                       "LinkerStretching.cu");
    CUDAcommon::handleerror(cudaMemcpy(gpu_kstr, kstr, numInteractions * sizeof(double), cudaMemcpyHostToDevice),
                            "cuda data transfer", "LinkerStretching.cu");

    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_eql, numInteractions * sizeof(double)),"cuda data transfer",
                            "LinkerStretching.cu");
    CUDAcommon::handleerror(cudaMemcpy(gpu_eql, eql, numInteractions * sizeof(double), cudaMemcpyHostToDevice),
                            "cuda data transfer", "LinkerStretching.cu");
    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_pos1, numInteractions * sizeof(double)),"cuda data transfer",
                            "LinkerStretching.cu");
    CUDAcommon::handleerror(cudaMemcpy(gpu_pos1, pos1, numInteractions * sizeof(double), cudaMemcpyHostToDevice),
                            "cuda data transfer", "LinkerStretching.cu");
    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_pos2, numInteractions * sizeof(double)),"cuda data transfer",
                            "LinkerStretching.cu");
    CUDAcommon::handleerror(cudaMemcpy(gpu_pos2, pos2, numInteractions * sizeof(double), cudaMemcpyHostToDevice),
                            "cuda data transfer", "LinkerStretching.cu");
    vector<int> params;
    params.push_back(int(n));
    params.push_back(numInteractions);

    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_params, 2 * sizeof(int)),"cuda data transfer",
                            "LinkerStretching.cu");
    CUDAcommon::handleerror(cudaMemcpy(gpu_params, params.data(), 2 * sizeof(int), cudaMemcpyHostToDevice),
                            "cuda data transfer", "LinkerStretching.cu");
    nvtxRangePop();
#endif
}

template<class LStretchingInteractionType>
void LinkerStretching<LStretchingInteractionType>::deallocate() {

    delete beadSet;
    delete kstr;
    delete eql;
    delete pos1;
    delete pos2;
#ifdef CUDAACCL
    _FFType.deallocate();
    CUDAcommon::handleerror(cudaFree(gpu_beadSet),"cudaFree", "LinkerStretching.cu");
    CUDAcommon::handleerror(cudaFree(gpu_kstr),"cudaFree", "LinkerStretching.cu");
    CUDAcommon::handleerror(cudaFree(gpu_pos1),"cudaFree", "LinkerStretching.cu");
    CUDAcommon::handleerror(cudaFree(gpu_pos2),"cudaFree", "LinkerStretching.cu");
    CUDAcommon::handleerror(cudaFree(gpu_eql),"cudaFree", "LinkerStretching.cu");
    CUDAcommon::handleerror(cudaFree(gpu_params),"cudaFree", "LinkerStretching.cu");
#endif
}


template <class LStretchingInteractionType>
double LinkerStretching<LStretchingInteractionType>::computeEnergy(double* coord, double *f, double d){

    double U_i[1], U_ii;
    double* gU_i;
    U_ii = NULL;
#ifdef CUDAACCL
    //has to be changed to accomodate aux force
    double * gpu_coord=CUDAcommon::getCUDAvars().gpu_coord;
    double * gpu_force=CUDAcommon::getCUDAvars().gpu_force;
    double * gpu_d = CUDAcommon::getCUDAvars().gpu_lambda;
    nvtxRangePushA("CCEL");

//    if(d == 0.0){
//        gU_i=_FFType.energy(gpu_coord, gpu_force, gpu_beadSet, gpu_kstr, gpu_eql, gpu_pos1, gpu_pos2, gpu_params);
//    }
//    else{
        gU_i=_FFType.energy(gpu_coord, gpu_force, gpu_beadSet, gpu_kstr, gpu_eql, gpu_pos1, gpu_pos2, gpu_d,
                            gpu_params);
//    }
    nvtxRangePop();
#else
    nvtxRangePushA("SCEL");

    if (d == 0.0)
        U_ii = _FFType.energy(coord, f, beadSet, kstr, eql, pos1, pos2);
    else
        U_ii = _FFType.energy(coord, f, beadSet, kstr, eql, pos1, pos2, d);
//    std::cout<<"================="<<endl;
    nvtxRangePop();
#endif

    return U_ii;
}

template <class LStretchingInteractionType>
void LinkerStretching<LStretchingInteractionType>::computeForces(double *coord, double *f) {
#ifdef CUDAACCL
    //has to be changed to accomodate aux force
    double * gpu_coord=CUDAcommon::getCUDAvars().gpu_coord;

    double * gpu_force;

//    //TODO remove this later need not copy forces back to CPU.
//    CUDAcommon::handleerror(cudaMemcpy(F_c, gpu_force, 3 * Bead::getBeads().size() *sizeof(double),
//                                       cudaMemcpyDeviceToHost));
//    cout.precision(dbl::max_digits10);
//    for(int iter=0;iter<Bead::getBeads().size();iter++) {
//        std::cout << "C " << F_c[3 * iter] << " " << F_c[3 * iter + 1] << " " << F_c[3 * iter + 2] <<" ";
//        std::cout << "V "<<f[3 * iter] << " " << f[3 * iter + 1] << " " << f[3 * iter + 2] << endl;
//    }
//    std::cout<<"check ends "<<blocksnthreads[0]<<" "<<blocksnthreads[1]<<endl;

    if(cross_checkclass::Aux){
        nvtxRangePushA("CCFL");

        gpu_force=CUDAcommon::getCUDAvars().gpu_forceAux;
        _FFType.forces(gpu_coord, gpu_force, gpu_beadSet, gpu_kstr, gpu_eql, gpu_pos1, gpu_pos2, gpu_params);
        nvtxRangePop();
    }
    else {
        nvtxRangePushA("CCFL");

        gpu_force = CUDAcommon::getCUDAvars().gpu_force;
        _FFType.forces(gpu_coord, gpu_force, gpu_beadSet, gpu_kstr, gpu_eql, gpu_pos1, gpu_pos2, gpu_params);
        nvtxRangePop();
    }

    //TODO remove this later need not copy forces back to CPU.
//    CUDAcommon::handleerror(cudaMemcpy(F_i, gpu_force, 3 * Bead::getBeads().size() *sizeof(double),
//                                       cudaMemcpyDeviceToHost),"cuda data transfer", "LinkerStretching.cu");
#else
    nvtxRangePushA("SCFL");
    _FFType.forces(coord, f, beadSet, kstr, eql, pos1, pos2);
    nvtxRangePop();
#endif
}


///Temlate specializations
template double LinkerStretching<LinkerStretchingHarmonic>::computeEnergy(double *coord, double *f, double d);
template void LinkerStretching<LinkerStretchingHarmonic>::computeForces(double *coord, double *f);
template void LinkerStretching<LinkerStretchingHarmonic>::vectorize();
template void LinkerStretching<LinkerStretchingHarmonic>::deallocate();

