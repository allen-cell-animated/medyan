
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

#include "FilamentStretching.h"

#include "FilamentStretchingHarmonic.h"
#include "Bead.h"
#include "cross_check.h"
#include "nvToolsExt.h"

template <class FStretchingInteractionType>
void FilamentStretching<FStretchingInteractionType>::vectorize() {

    beadSet = new int[n * Cylinder::getCylinders().size()];
    kstr = new double[Cylinder::getCylinders().size()];
    eql = new double[Cylinder::getCylinders().size()];

    int i = 0;

    for (auto c: Cylinder::getCylinders()) {
        beadSet[n * i] = c->getFirstBead()->_dbIndex;
        beadSet[n * i + 1] = c->getSecondBead()->_dbIndex;

        kstr[i] = c->getMCylinder()->getStretchingConst();
        eql[i] = c->getMCylinder()->getEqLength();

        i++;
    }
    //CUDA
#ifdef CUDAACCL
    nvtxRangePushA("CVFF");
//    F_i = new double[3 * Bead::getBeads().size()];
    int numInteractions = Cylinder::getCylinders().size();
    _FFType.optimalblocksnthreads(numInteractions);

    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_beadSet, n * numInteractions * sizeof(int)),"cuda data "
            "transfer", " FilamentStretching.cu");
    CUDAcommon::handleerror(cudaMemcpy(gpu_beadSet, beadSet, n * numInteractions * sizeof(int),
                                       cudaMemcpyHostToDevice),"cuda data transfer", " FilamentStretching.cu");

    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_kstr, numInteractions * sizeof(double)),"cuda data transfer",
                            " FilamentStretching.cu");
    CUDAcommon::handleerror(cudaMemcpy(gpu_kstr, kstr, numInteractions * sizeof(double), cudaMemcpyHostToDevice),
                            "cuda data transfer", " FilamentStretching.cu");

    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_eql, numInteractions * sizeof(double)),"cuda data transfer",
                            " FilamentStretching.cu");
    CUDAcommon::handleerror(cudaMemcpy(gpu_eql, eql, numInteractions * sizeof(double), cudaMemcpyHostToDevice),
                            "cuda data transfer", " FilamentStretching.cu");

    vector<int> params;
    params.push_back(int(n));
    params.push_back(numInteractions);
    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_params, 2 * sizeof(int)),"cuda data transfer",
                            " FilamentStretching.cu");
    CUDAcommon::handleerror(cudaMemcpy(gpu_params, params.data(), 2 * sizeof(int), cudaMemcpyHostToDevice),
            "cuda data transfer", " FilamentStretching.cu");
    nvtxRangePop();
#endif

    //
}

template<class FStretchingInteractionType>
void FilamentStretching<FStretchingInteractionType>::deallocate() {

    delete beadSet;
    delete kstr;
    delete eql;
#ifdef CUDAACCL
    _FFType.deallocate();
    CUDAcommon::handleerror(cudaFree(gpu_beadSet));
    CUDAcommon::handleerror(cudaFree(gpu_kstr));
    CUDAcommon::handleerror(cudaFree(gpu_eql));
    CUDAcommon::handleerror(cudaFree(gpu_params));
#endif
}


template <class FStretchingInteractionType>
double FilamentStretching<FStretchingInteractionType>::computeEnergy(double* coord, double *f, double d){

    double U_i[1], U_ii;
    double* gU_i;
    U_ii = NULL;
#ifdef CUDAACCL
    //has to be changed to accomodate aux force
    double * gpu_coord=CUDAcommon::getCUDAvars().gpu_coord;
    double * gpu_force=CUDAcommon::getCUDAvars().gpu_force;
    double * gpu_d = CUDAcommon::getCUDAvars().gpu_lambda;
//    std::cout<<"Fil Stretching Forces"<<endl;
    nvtxRangePushA("CCEFS");

//    if(d == 0.0){
//        gU_i=_FFType.energy(gpu_coord, gpu_force, gpu_beadSet, gpu_kstr, gpu_eql, gpu_params);
//
//    }
//    else{
        gU_i=_FFType.energy(gpu_coord, gpu_force, gpu_beadSet, gpu_kstr, gpu_eql, gpu_d,
                            gpu_params);
//    }
    nvtxRangePop();
#else
    nvtxRangePushA("SCEFS");

    if (d == 0.0)
        U_ii = _FFType.energy(coord, f, beadSet, kstr, eql);
    else
        U_ii = _FFType.energy(coord, f, beadSet, kstr, eql, d);
    nvtxRangePop();
#endif
//    whoisCulprit();
    return U_ii;
}

template <class FStretchingInteractionType>
void FilamentStretching<FStretchingInteractionType>::computeForces(double *coord, double *f) {
#ifdef CUDAACCL
    //has to be changed to accomodate aux force
    double * gpu_coord=CUDAcommon::getCUDAvars().gpu_coord;

    double * gpu_force;
    if(cross_checkclass::Aux){
        nvtxRangePushA("CCFFS");

        gpu_force=CUDAcommon::getCUDAvars().gpu_forceAux;
        _FFType.forces(gpu_coord, gpu_force, gpu_beadSet, gpu_kstr, gpu_eql, gpu_params);
        nvtxRangePop();
    }
    else {
        nvtxRangePushA("CCFFS");
        gpu_force = CUDAcommon::getCUDAvars().gpu_force;
        _FFType.forces(gpu_coord, gpu_force, gpu_beadSet, gpu_kstr, gpu_eql, gpu_params);
        nvtxRangePop();
    }

#else
    nvtxRangePushA("SCFFS");
    _FFType.forces(coord, f, beadSet, kstr, eql);
    nvtxRangePop();
#endif

}

//template <class FStretchingInteractionType>
//void FilamentStretching<FStretchingInteractionType>::whoisCulprit() {
//    cudaDeviceSynchronize();
//    std::cout<<cudaGetLastError()<<endl;
//    cout<<cudaGetErrorString(cudaGetLastError())<<endl;
//    if( cudaGetLastError() == cudaErrorAssert){
//        cudaDeviceSynchronize();
//        auto culpritinteraction = CUDAcommon::getCUDAvars().culpritinteraction;
//        std::cout<<strcmp(culpritinteraction,"Filament Stretching")<<endl;
//        if(strcmp(culpritinteraction, "Filament Stretching")==0){
//            CUDAcommon::printculprit("FilamentStretching","FilamentStretchingHarmonic");
//            Filament* f;
//            f = (Filament*)(Cylinder::getCylinders()[CUDAcommon::getCUDAvars().culpritID[0]]->getParent());
//            cout<<"Printing culprit Filament information."<<endl;
//            f->printSelf();
//            exit(EXIT_FAILURE);
//        }
//    }
//}
///Temlate specializations
template double FilamentStretching<FilamentStretchingHarmonic>::computeEnergy(double *coord, double *f, double d);
template void FilamentStretching<FilamentStretchingHarmonic>::computeForces(double *coord, double *f);
template void FilamentStretching<FilamentStretchingHarmonic>::vectorize();
template void FilamentStretching<FilamentStretchingHarmonic>::deallocate();
//template void FilamentStretching<FilamentStretchingHarmonic>::whoisCulprit();
