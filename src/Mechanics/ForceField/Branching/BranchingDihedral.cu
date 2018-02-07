
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

#include "BranchingDihedral.h"

#include "BranchingDihedralCosine.h"

#include "BranchingPoint.h"
#include "Cylinder.h"
#include "Bead.h"
#include "cross_check.h"
#include "nvToolsExt.h"

template <class BDihedralInteractionType>
void BranchingDihedral<BDihedralInteractionType>::vectorize() {

    beadSet = new int[n * BranchingPoint::getBranchingPoints().size()];
    kdih = new double[BranchingPoint::getBranchingPoints().size()];
    pos = new double[BranchingPoint::getBranchingPoints().size()];

    int i = 0;

    for (auto b: BranchingPoint::getBranchingPoints()) {

        beadSet[n * i] = b->getFirstCylinder()->getFirstBead()->_dbIndex;
        beadSet[n * i + 1] = b->getFirstCylinder()->getSecondBead()->_dbIndex;
        beadSet[n * i + 2] = b->getSecondCylinder()->getFirstBead()->_dbIndex;
        beadSet[n * i + 3] = b->getSecondCylinder()->getSecondBead()->_dbIndex;

        kdih[i] = b->getMBranchingPoint()->getDihedralConstant();
        pos[i] = b->getPosition();

        i++;
    }
    //CUDA
#ifdef CUDAACCL
//    cudaEvent_t start, stop;
//    CUDAcommon::handleerror(cudaEventCreate( &start));
//    CUDAcommon::handleerror(cudaEventCreate( &stop));
//    CUDAcommon::handleerror(cudaEventRecord( start, 0));
    nvtxRangePushA("CVFF");

    int numInteractions =BranchingPoint::getBranchingPoints().size();
    _FFType.optimalblocksnthreads(numInteractions);
    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_beadSet, n * numInteractions * sizeof(int)));
    CUDAcommon::handleerror(cudaMemcpy(gpu_beadSet, beadSet, n * numInteractions * sizeof(int),
                                       cudaMemcpyHostToDevice));

    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_kdih, numInteractions * sizeof(double)));
    CUDAcommon::handleerror(cudaMemcpy(gpu_kdih, kdih, numInteractions * sizeof(double), cudaMemcpyHostToDevice));

    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_pos, numInteractions * sizeof(double)));
    CUDAcommon::handleerror(cudaMemcpy(gpu_pos, pos, numInteractions * sizeof(double), cudaMemcpyHostToDevice));

    vector<int> params;
    params.push_back(int(n));
    params.push_back(numInteractions);

    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_params, 2 * sizeof(int)));
    CUDAcommon::handleerror(cudaMemcpy(gpu_params, params.data(), 2 * sizeof(int), cudaMemcpyHostToDevice));

    nvtxRangePop();

#endif
}

template<class BDihedralInteractionType>
void BranchingDihedral<BDihedralInteractionType>::deallocate() {

    delete beadSet;
    delete kdih;
    delete pos;
#ifdef CUDAACCL
    _FFType.deallocate();
    CUDAcommon::handleerror(cudaFree(gpu_beadSet));
    CUDAcommon::handleerror(cudaFree(gpu_kdih));
    CUDAcommon::handleerror(cudaFree(gpu_pos));
    CUDAcommon::handleerror(cudaFree(gpu_params));
#endif
}


template <class BDihedralInteractionType>
double BranchingDihedral<BDihedralInteractionType>::computeEnergy(double *coord, double *f, double d) {

    double U_i[1], U_ii;
    double* gU_i;
    U_ii = NULL;
#ifdef CUDAACCL
    //has to be changed to accomodate aux force
    double * gpu_coord=CUDAcommon::getCUDAvars().gpu_coord;
    double * gpu_force=CUDAcommon::getCUDAvars().gpu_force;
    double * gpu_d = CUDAcommon::getCUDAvars().gpu_lambda;
    nvtxRangePushA("CCEBD");

    if(d == 0.0){
        gU_i=_FFType.energy(gpu_coord, gpu_force, gpu_beadSet, gpu_kdih, gpu_pos, gpu_params);

    }
    else{
        gU_i=_FFType.energy(gpu_coord, gpu_force, gpu_beadSet, gpu_kdih, gpu_pos, gpu_d,
                            gpu_params);
    }
    nvtxRangePop();
#endif
    nvtxRangePushA("SCEBD");

    if (d == 0.0)
        U_ii = _FFType.energy(coord, f, beadSet, kdih, pos);
    else
        U_ii = _FFType.energy(coord, f, beadSet, kdih, pos, d);
    nvtxRangePop();
    if(gU_i!=NULL) {

        CUDAcommon::handleerror(cudaMemcpy(U_i, gU_i, sizeof(double),
                                           cudaMemcpyDeviceToHost));
    }
    else
        U_i[0] = 0.0;
    if(fabs(U_ii)>1000000.0) {
        if (fabs((U_ii - U_i[0]) / U_ii) > 0.0001){
            std::cout<<endl;
            std::cout << "CUDA BDE " << U_i[0] << endl;
            std::cout << "Vectorized BDE " << U_ii << endl;
            std::cout << "Precision match error" << fabs(U_ii - U_i[0]) << endl;
        }
    }
    else {
        if (fabs(U_ii - U_i[0]) > 1.0 / 100000000.0){
            std::cout<<endl;
            std::cout << "CUDA BDE " << U_i << endl;
            std::cout << "Vectorized BDE " << U_ii << endl;
            std::cout << "Precision match " << fabs(U_ii - U_i[0]) << endl;
        }
    }

    return U_ii;

}

template <class BDihedralInteractionType>
void BranchingDihedral<BDihedralInteractionType>::computeForces(double *coord, double *f) {
#ifdef CUDAACCL
    //has to be changed to accomodate aux force
    double * gpu_coord=CUDAcommon::getCUDAvars().gpu_coord;

    double * gpu_force;

    if(cross_checkclass::Aux){
        nvtxRangePushA("CCFBD");

        gpu_force=CUDAcommon::getCUDAvars().gpu_forceAux;
        _FFType.forces(gpu_coord, gpu_force, gpu_beadSet, gpu_kdih, gpu_pos, gpu_params);
        nvtxRangePop();
    }
    else {
        nvtxRangePushA("CCFBD");

        gpu_force = CUDAcommon::getCUDAvars().gpu_force;
        _FFType.forces(gpu_coord, gpu_force, gpu_beadSet, gpu_kdih, gpu_pos, gpu_params);
        nvtxRangePop();
    }

    //TODO remove this later need not copy forces back to CPU.
    CUDAcommon::handleerror(cudaMemcpy(F_i, gpu_force, 3 * Bead::getBeads().size() *sizeof(double),
                                       cudaMemcpyDeviceToHost));
#endif
    nvtxRangePushA("SCFBD");

    _FFType.forces(coord, f, beadSet, kdih, pos);
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
            std::cout<<"BD Forces"<<endl;
            std::cout << "CUDA       " << F_i[3 * iter] << " " << F_i[3 * iter + 1] << " " << F_i[3 * iter + 2] << endl;
            std::cout << "Vectorized " << f[3 * iter] << " " << f[3 * iter + 1] << " " << f[3 * iter + 2] << endl;
            std::cout<<"Precision match "<<fabs(F_i[3 * iter] - f[3 * iter])<<" "<<fabs(F_i[3 * iter + 1] - f[3 *
                                                                                                              iter + 1])<<" "<<fabs(F_i[3 * iter + 2] - f[3 * iter + 2])<<endl;
//        exit(EXIT_FAILURE);
        }
    }
#endif
}

///Template specializations
template double BranchingDihedral<BranchingDihedralCosine>::computeEnergy(double *coord, double *f, double d);
template void BranchingDihedral<BranchingDihedralCosine>::computeForces(double *coord, double *f);
template void BranchingDihedral<BranchingDihedralCosine>::vectorize();
template void BranchingDihedral<BranchingDihedralCosine>::deallocate();
