
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

#include <src/Mechanics/Minimizer/CGMethod.h>
#include "BranchingStretching.h"

#include "BranchingStretchingHarmonic.h"

#include "BranchingPoint.h"
#include "Cylinder.h"
#include "Bead.h"
#include "cross_check.h"
#include "nvToolsExt.h"

template <class BStretchingInteractionType>
void BranchingStretching<BStretchingInteractionType>::vectorize() {

    beadSet = new int[n * BranchingPoint::getBranchingPoints().size()];
    kstr = new double[BranchingPoint::getBranchingPoints().size()];
    eql = new double[BranchingPoint::getBranchingPoints().size()];
    pos = new double[BranchingPoint::getBranchingPoints().size()];


    int i = 0;

    for (auto b: BranchingPoint::getBranchingPoints()) {

        beadSet[n * i] = b->getFirstCylinder()->getFirstBead()->_dbIndex;
        beadSet[n * i + 1] = b->getFirstCylinder()->getSecondBead()->_dbIndex;
        beadSet[n * i + 2] = b->getSecondCylinder()->getFirstBead()->_dbIndex;

        kstr[i] = b->getMBranchingPoint()->getStretchingConstant();
        eql[i] = b->getMBranchingPoint()->getEqLength();
        pos[i] = b->getPosition();

        i++;
    }
    //CUDA
#ifdef CUDAACCL
    F_i = new double[CGMethod::N];
    nvtxRangePushA("CVFF");

    int numInteractions = BranchingPoint::getBranchingPoints().size();
    _FFType.optimalblocksnthreads(numInteractions);

    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_beadSet, n * numInteractions * sizeof(int)));
    CUDAcommon::handleerror(cudaMemcpy(gpu_beadSet, beadSet, n * numInteractions * sizeof(int),
                                       cudaMemcpyHostToDevice));

    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_kstr, numInteractions * sizeof(double)));
    CUDAcommon::handleerror(cudaMemcpy(gpu_kstr, kstr, numInteractions * sizeof(double), cudaMemcpyHostToDevice));

    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_eql, numInteractions * sizeof(double)));
    CUDAcommon::handleerror(cudaMemcpy(gpu_eql, eql, numInteractions * sizeof(double), cudaMemcpyHostToDevice));

    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_pos, numInteractions * sizeof(double)));
    CUDAcommon::handleerror(cudaMemcpy(gpu_pos, pos, numInteractions * sizeof(double), cudaMemcpyHostToDevice));

    vector<int> params;
    params.push_back(int(n));
    params.push_back(numInteractions);
    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_params, 2 * sizeof(int)));
    CUDAcommon::handleerror(cudaMemcpy(gpu_params, params.data(), 2 * sizeof(int), cudaMemcpyHostToDevice));
    nvtxRangePop();
#endif

    //
}

template<class BStretchingInteractionType>
void BranchingStretching<BStretchingInteractionType>::deallocate() {

    delete beadSet;
    delete kstr;
    delete eql;
    delete pos;
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
double BranchingStretching<BStretchingInteractionType>::computeEnergy(double *coord, double *f, double d) {


    double U_i[1], U_ii;
    double* gU_i;
    U_ii = NULL;
#ifdef CUDAACCL
    //has to be changed to accomodate aux force
    double * gpu_coord=CUDAcommon::getCUDAvars().gpu_coord;
    double * gpu_force=CUDAcommon::getCUDAvars().gpu_force;
    double * gpu_d = CUDAcommon::getCUDAvars().gpu_lambda;
    nvtxRangePushA("CCEBS");

    if(d == 0.0){
        gU_i=_FFType.energy(gpu_coord, gpu_force, gpu_beadSet, gpu_kstr, gpu_eql, gpu_pos, gpu_params);

    }
    else{
        gU_i=_FFType.energy(gpu_coord, gpu_force, gpu_beadSet, gpu_kstr, gpu_eql, gpu_pos, gpu_d,
                            gpu_params);
    }
    nvtxRangePop();
#endif
    nvtxRangePushA("SCEBS");
    if (d == 0.0)
        U_ii = _FFType.energy(coord, f, beadSet, kstr, eql, pos);
    else
        U_ii = _FFType.energy(coord, f, beadSet, kstr, eql, pos, d);
    nvtxRangePop();
    if(gU_i!=NULL) {

        CUDAcommon::handleerror(cudaMemcpy(U_i, gU_i, sizeof(double), cudaMemcpyDeviceToHost),"computeEnergy",
                                "BranchingStretching.cu");
    }
    else
        U_i[0] = 0.0;
    if(fabs(U_ii)>1000000.0) {
        if (fabs((U_ii - U_i[0]) / U_ii) > 0.0001){
            std::cout<<endl;
            std::cout << "CUDA BSE " << U_i[0] << endl;
            std::cout << "Vectorized BSE " << U_ii << endl;
            std::cout << "Precision match error" << fabs(U_ii - U_i[0]) << endl;
        }
    }
    else {
        if (fabs(U_ii - U_i[0]) > 1.0 / 100000000.0){
            std::cout<<endl;
            std::cout << "CUDA BSE " << U_i[0] << endl;
            std::cout << "Vectorized BSE " << U_ii << endl;
            std::cout << "Precision match " << fabs(U_ii - U_i[0]) << endl;
//        exit(EXIT_FAILURE);
        }
    }

    return U_ii;
}

template <class BStretchingInteractionType>
void BranchingStretching<BStretchingInteractionType>::computeForces(double *coord, double *f) {
#ifdef CUDAACCL
    //has to be changed to accomodate aux force
    double * gpu_coord=CUDAcommon::getCUDAvars().gpu_coord;

    double * gpu_force;


    if(cross_checkclass::Aux){
        nvtxRangePushA("CCFBS");

        gpu_force=CUDAcommon::getCUDAvars().gpu_forceAux;
        _FFType.forces(gpu_coord, gpu_force, gpu_beadSet, gpu_kstr, gpu_eql, gpu_pos, gpu_params);
        nvtxRangePop();
    }
    else {
        nvtxRangePushA("CCFBS");

        gpu_force = CUDAcommon::getCUDAvars().gpu_force;
        _FFType.forces(gpu_coord, gpu_force, gpu_beadSet, gpu_kstr, gpu_eql, gpu_pos, gpu_params);
        nvtxRangePop();
    }
    //TODO remove this later need not copy forces back to CPU.
    CUDAcommon::handleerror(cudaMemcpy(F_i, gpu_force, CGMethod::N*sizeof(double),
                                       cudaMemcpyDeviceToHost));
#endif
    nvtxRangePushA("SCFBS");

    _FFType.forces(coord, f, beadSet, kstr, eql, pos);
    nvtxRangePop();
#ifdef CUDAACCL

    bool state = false;
    for(auto iter=0;iter<Bead::getBeads().size();iter++) {
        if (fabs(F_i[3 * iter] - f[3 * iter]) <=1.0/100000000.0 && fabs(F_i[3 * iter + 1] - f[3 * iter + 1])
                                                                   <=1.0/100000000.0 && fabs(F_i[3 * iter + 2] - f[3 * iter + 2]) <=1.0/100000000.0)
        {state = true;}
        else {
            state = false;
            std::cout<<endl;
            std::cout<<"BS Forces"<<endl;
            std::cout<<"Precision match "<<fabs(F_i[3 * iter] - f[3 * iter])<<" "<<fabs(F_i[3 * iter + 1] - f[3 *
                                                                                                              iter + 1])<<" "<<fabs(F_i[3 * iter + 2] - f[3 * iter + 2])<<endl;
            std::cout << "CUDA       " << F_i[3 * iter] << " " << F_i[3 * iter + 1] << " " << F_i[3 * iter + 2] << endl;
            std::cout << "Vectorized " << f[3 * iter] << " " << f[3 * iter + 1] << " " << f[3 * iter + 2] << endl;
//        exit(EXIT_FAILURE);
        }
    }
#endif
}



///Template specializations
template double
BranchingStretching<BranchingStretchingHarmonic>::computeEnergy(double *coord, double *f, double d);
template void BranchingStretching<BranchingStretchingHarmonic>::computeForces(double *coord, double *f);
template void BranchingStretching<BranchingStretchingHarmonic>::vectorize();
template void BranchingStretching<BranchingStretchingHarmonic>::deallocate();
