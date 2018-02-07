
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

#include "BoundaryCylinderRepulsionExp.h"
#include "BoundaryCylinderRepulsion.h"

#include "BoundaryElement.h"
#include "Bead.h"
//TODO added for temporary CUDA force.
#include "CGMethod.h"
#include "BoundaryCylinderRepulsionExpCUDA.h"
#include "CUDAcommon.h"
#include "cross_check.h"
double BoundaryCylinderRepulsionExp::energy(double *coord, double *f, int *beadSet,
                                            double *krep, double *slen, int *nneighbors) {

    int nb, nc;
    double *coord1, R, r, U_i;
    double U = 0;
    auto Cumnc=0;
    auto beList = BoundaryElement::getBoundaryElements();
    nb = beList.size();

    for (int ib = 0; ib < nb; ib++) {

        auto be = beList[ib];
        nc = nneighbors[ib];

        for (int ic = 0; ic < nc; ic++) {

            coord1 = &coord[3 * beadSet[Cumnc + ic]];
            r = be->distance(coord1);

            R = -r / slen[Cumnc + ic];
            U_i = krep[Cumnc + ic] * exp(R);

//            std::cout<<r<<" "<<U_i<<endl;
            if (fabs(U_i) == numeric_limits<double>::infinity()
                || U_i != U_i || U_i < -1.0) {

                //set culprit and return
                BoundaryInteractions::_boundaryElementCulprit = be;
                ///TODO
                //BoundaryInteractions::_otherCulprit;

                return -1;
            }
            U += U_i;
        }
        Cumnc += nc;
    }
    return U;
}

double BoundaryCylinderRepulsionExp::energy(double *coord, double *f, int *beadSet,
                                            double *krep, double *slen, int *nneighbors, double d) {

    int nb, nc;
    double *coord1, *force1, R, r, U_i;
    double U = 0;
    long Cumnc=0;
    auto beList = BoundaryElement::getBoundaryElements();
    nb = beList.size();

    for (int ib = 0; ib < nb; ib++) {

        auto be = beList[ib];
        nc = nneighbors[ib];

        for(int ic = 0; ic < nc; ic++) {

            coord1 = &coord[3 * beadSet[Cumnc + ic]];
            force1 = &f[3 * beadSet[Cumnc + ic]];

            r = be->stretchedDistance(coord1, force1, d);

            R = -r / slen[Cumnc + ic];

//            std::cout<<r<<" "<<krep[Cumnc+ic]<<endl;
            U_i = krep[Cumnc + ic] * exp(R);


            if(fabs(U_i) == numeric_limits<double>::infinity()
               || U_i != U_i || U_i < -1.0) {

                //set culprit and return
                BoundaryInteractions::_boundaryElementCulprit = be;
                ///TODO
                //BoundaryInteractions::_otherCulprit;

                return -1;
            }
            U += U_i;
        }
        Cumnc+=nc;
    }
    return U;
}



void BoundaryCylinderRepulsionExp::forces(double *coord, double *f, int *beadSet,
                                          double *krep, double *slen, int *nneighbors) {
    int nb, nc;
    double *coord1, *force1, R, r, f0;
    double *F_i;
    double *forcecopy;
    std::cout<<"bdry f init vec"<<endl;
    forcecopy = new double[CGMethod::N];
    for(auto iter=0;iter<CGMethod::N;iter++)
        forcecopy[iter]=0.0;

    auto beList = BoundaryElement::getBoundaryElements();
    nb = beList.size();
    auto Cumnc=0;
//     for (int ib = 0; ib < nb; ib++) {
//         auto be = beList[ib];
//         nc = nneighbors[ib];
//         std::cout<<ib<<" "<<nc<<endl;
//          for(int ic = 0; ic < nc; ic++) {
//              std::cout<<beadSet[ Cumnc + ic]<<" ";
//          }
//         std::cout<<endl;
//         Cumnc+=nc;
//
//     }
//    for(auto i=0;i<CGMethod::N;i++)
//        std::cout<<f[i]<<" ";
//    std::cout<<endl;

    for (int ib = 0; ib < nb; ib++) {

        auto be = beList[ib];
        nc = nneighbors[ib];
        for(int ic = 0; ic < nc; ic++) {
            coord1 = &coord[3 * beadSet[ Cumnc + ic]];
            force1 = &f[3 * beadSet[ Cumnc + ic]];
            r = be->distance(coord1);
            auto norm = be->normal(coord1);

            R = -r / slen[Cumnc + ic];
            f0 = krep[Cumnc + ic] * exp(R);

            force1[0] += f0 *norm[0];
            force1[1] += f0 *norm[1];
            force1[2] += f0 *norm[2];

            forcecopy[3 * beadSet[ Cumnc + ic]] += f0 *norm[0];
            forcecopy[3 * beadSet[ Cumnc + ic] + 1] += f0 *norm[1];
            forcecopy[3 * beadSet[ Cumnc + ic] + 2] += f0 *norm[2];

            //            std::cout<<beadSet[Cumnc+ic]<<"
            // "<<norm[0]<<" "<<norm[1]<<"
            // "<<norm[2]<<" "<<f0<<endl;
//                        std::cout<<beadSet[ Cumnc + ic]<<" "<<force1[0]<<" "<<force1[1]<<" "<<force1[2]<<" "<<slen[Cumnc+ic]<<" "
//                        <<krep[Cumnc+ic]<<" "<<f0<<endl;

        }
        Cumnc+=nc;
    }

    double *gfcopy;
    double *gpu_force;
    vector<int> blocksnthreads;
    int *gpu_nint; int nint[1]; nint[0]=CGMethod::N/3;
    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_nint, sizeof(int)));
    CUDAcommon::handleerror(cudaMemcpy(gpu_nint, nint, sizeof(int), cudaMemcpyHostToDevice));
    blocksnthreads.push_back(CGMethod::N/(3*THREADSPERBLOCK) + 1);
    if(blocksnthreads[0]==1) blocksnthreads.push_back(CGMethod::N/3);
    else blocksnthreads.push_back(THREADSPERBLOCK);
//    std::cout<<"Cpy Bdry Frc Number of Blocks: "<<blocksnthreads[0]<<endl;
//    std::cout<<"Threads per block: "<<blocksnthreads[1]<<endl;
    CUDAcommon::handleerror(cudaMalloc((void **) &gfcopy, CGMethod::N * sizeof(double)));
    CUDAcommon::handleerror(cudaMemcpy(gfcopy, forcecopy, CGMethod::N * sizeof(double), cudaMemcpyHostToDevice));
    if(cross_checkclass::Aux) {
        gpu_force = CUDAcommon::getCUDAvars().gpu_forceAux;
    }
    else
        gpu_force = CUDAcommon::getCUDAvars().gpu_force;

    //TODO remove this later need not copy forces back to CPU.
//    CUDAcommon::handleerror(cudaMemcpy(F_i, gpu_force, 3 * Bead::getBeads().size() *sizeof(double),
//                                       cudaMemcpyDeviceToHost));
//    for(auto i=0;i<CGMethod::N;i++)
//        std::cout<<F_i[i]<<" ";
//    std::cout<<endl;
    cudaStream_t  stream;
    CUDAcommon::handleerror(cudaStreamCreate(&stream));
    BoundaryCylinderRepulsionadd << < blocksnthreads[0], blocksnthreads[1],0, stream>> >(gpu_force, gfcopy, gpu_nint);
    CUDAcommon::handleerror(cudaGetLastError(),"BoundaryCylinderRepulsionadd", "BoundaryCylinderRepulsionExp.cu");
    CUDAcommon::handleerror(cudaStreamSynchronize(stream),"BoundaryCylinderRepulsionadd", "BoundaryCylinderRepulsionExp.cu");
    CUDAcommon::handleerror(cudaStreamDestroy(stream),"BoundaryCylinderRepulsionadd", "BoundaryCylinderRepulsionExp.cu");
//    CUDAcommon::handleerror( cudaPeekAtLastError() );
//    CUDAcommon::handleerror(cudaDeviceSynchronize());
    std::cout<<"bdry f cpy 2 hostt"<<endl;
    F_i = new double[3*Bead::getBeads().size()];
    CUDAcommon::handleerror(cudaMemcpy(F_i, gpu_force, 3 * Bead::getBeads().size() *sizeof(double),
                                       cudaMemcpyDeviceToHost));
    CUDAcommon::handleerror(cudaFree(gfcopy));
#ifdef CUDAACCL
    cout.precision(dbl::max_digits10);
//    std::cout<<"B forces"<<endl;
//    for(int iter=0;iter<Bead::getBeads().size();iter++) {
//        std::cout << F_i[3 * iter] << " " << F_i[3 * iter + 1] << " " << F_i[3 * iter + 2] <<" ";
//        std::cout <<f[3 * iter] << " " << f[3 * iter + 1] << " " << f[3 * iter + 2] << endl;
//    }

    bool state = false;
    for(int iter=0;iter<Bead::getBeads().size();iter++) {
        if (fabs(F_i[3 * iter] - f[3 * iter]) <=1.0/100000000.0 && fabs(F_i[3 * iter + 1] - f[3 * iter + 1])
        <=1.0/100000000.0 && fabs(F_i[3 * iter + 2] - f[3 * iter + 2]) <=1.0/100000000.0)
        {state = true;}
        else {
            state = false;
            std::cout<<endl;
            std::cout<<"Precision match "<<fabs(F_i[3 * iter] - f[3 * iter])<<" "<<fabs(F_i[3 * iter + 1] - f[3 *
                                                                                                              iter + 1])<<" "<<fabs(F_i[3 * iter + 2] - f[3 * iter + 2])<<endl;
            std::cout << "CUDA       " << F_i[3 * iter] << " " << F_i[3 * iter + 1] << " " << F_i[3 * iter + 2] << endl;
            std::cout << "Vectorized " << f[3 * iter] << " " << f[3 * iter + 1] << " " << f[3 * iter + 2] << endl;

//        exit(EXIT_FAILURE);
        }
    }
//    if(state)
//    std::cout<<"F M+V+B YES"<<endl;
    delete [] F_i;
    delete [] forcecopy;
#endif
}

double BoundaryCylinderRepulsionExp::loadForces(double r, double kRep, double screenLength) {

    double R = -r/screenLength;
    return kRep * exp(R)/screenLength;

}