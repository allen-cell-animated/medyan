
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

#include <cmath>

#include "BranchingBendingCosine.h"
#include "BranchingBending.h"
#include "BranchingBendingCosineCUDA.h"

#include "BranchingPoint.h"
#include "Bead.h"

#include "MathFunctions.h"
#ifdef CUDAACCL
#include "nvToolsExt.h"
#endif

using namespace mathfunc;
#ifdef CUDAACCL
void BranchingBendingCosine::deallocate(){
    if(!(CUDAcommon::getCUDAvars().conservestreams))
        CUDAcommon::handleerror(cudaStreamDestroy(stream));
    CUDAcommon::handleerror(cudaFree(gU_i));
    CUDAcommon::handleerror(cudaFree(gU_sum));
    CUDAcommon::handleerror(cudaFree(gFF));
    CUDAcommon::handleerror(cudaFree(ginteraction));
}
void BranchingBendingCosine::optimalblocksnthreads( int nint){
    //CUDA stream create
    if(stream == NULL || !(CUDAcommon::getCUDAvars().conservestreams))
        CUDAcommon::handleerror(cudaStreamCreate(&stream));
    blocksnthreadse.clear();
    blocksnthreadsez.clear();
    blocksnthreadsf.clear();
    int blockSize;   // The launch configurator returned block size
    int minGridSize; // The minimum grid size needed to achieve the
    // maximum occupancy for a full device launch
    if(nint>0) {
        cudaOccupancyMaxPotentialBlockSizeVariableSMem(&minGridSize, &blockSize,
                                                       BranchingBendingCosineenergy, blockToSmem, 0);
        blocksnthreadse.push_back((nint + blockSize - 1) / blockSize);
        blocksnthreadse.push_back(blockSize);
        blockSize = 0;

        cudaOccupancyMaxPotentialBlockSizeVariableSMem(&minGridSize, &blockSize,
                                                       BranchingBendingCosineenergyz, blockToSmemez, 0);
        blocksnthreadsez.push_back((nint + blockSize - 1) / blockSize);
        blocksnthreadsez.push_back(blockSize);
        blockSize = 0;

        cudaOccupancyMaxPotentialBlockSizeVariableSMem(&minGridSize, &blockSize,
                                                       BranchingBendingCosineforces, blockToSmem, 0);
        blocksnthreadsf.push_back((nint + blockSize - 1) / blockSize);
        blocksnthreadsf.push_back(blockSize);

//        CUDAcommon::handleerror(cudaMalloc((void **) &gU_i, nint*sizeof(double)));
//        CUDAcommon::handleerror(cudaMalloc((void **) &gU_sum, sizeof(double)));
        //get addition vars
        bntaddvec2.clear();
        bntaddvec2 = getaddred2bnt(nint);
        CUDAcommon::handleerror(cudaMalloc((void **) &gU_i, bntaddvec2.at(0)*sizeof(double)));
        CUDAcommon::handleerror(cudaMemset(gU_i, 0, bntaddvec2.at(0) * sizeof(double)));
//        CUDAcommon::handleerror(cudaMalloc((void **) &gU_i, nint*sizeof(double)));
        CUDAcommon::handleerror(cudaMalloc((void **) &gU_sum, sizeof(double)));

        char a[] = "BranchingFF";
        char b[] = "Branching Bending Cosine";
        CUDAcommon::handleerror(cudaMalloc((void **) &gFF, 100 * sizeof(char)));
        CUDAcommon::handleerror(cudaMalloc((void **) &ginteraction, 100 * sizeof(char)));
        CUDAcommon::handleerror(cudaMemcpy(gFF, a, 100 * sizeof(char), cudaMemcpyHostToDevice));
        CUDAcommon::handleerror(cudaMemcpy(ginteraction, b, 100 * sizeof(char), cudaMemcpyHostToDevice));
    }
    else{
        blocksnthreadse.push_back(0);
        blocksnthreadse.push_back(0);
        blocksnthreadsez.push_back(0);
        blocksnthreadsez.push_back(0);
        blocksnthreadsf.push_back(0);
        blocksnthreadsf.push_back(0);
    }

}
double* BranchingBendingCosine::energy(double *coord, double *f, int *beadSet,
                                       double *kbend, double *eqt, int *params) {
//    if(blocksnthreadse[1]>0) {
//        BranchingBendingCosineenergy<<<blocksnthreadse[0], blocksnthreadse[1], (12 * blocksnthreadse[1]) * sizeof
//                (double), stream>>> (coord, f, beadSet, kbend, eqt, params, gU_i, CUDAcommon::getCUDAvars().gculpritID,
//                CUDAcommon::getCUDAvars().gculpritFF,
//                CUDAcommon::getCUDAvars().gculpritinteraction, gFF, ginteraction);
//        auto cvars = CUDAcommon::getCUDAvars();
//        cvars.streamvec.push_back(&stream);
//        CUDAcommon::cudavars = cvars;
//        CUDAcommon::handleerror( cudaGetLastError() ,"BranchingBendingCosineenergy", "BranchingBendingCosine.cu");
//        double* gpu_Utot = CUDAcommon::getCUDAvars().gpu_energy;
//        addvector<<<1,1,0,stream>>>(gU_i,params, gU_sum, gpu_Utot);
//        CUDAcommon::handleerror( cudaGetLastError(),"BranchingBendingCosineenergy", "BranchingBendingCosine.cu");
//        return gU_sum;}
//    else
//        return NULL;
}


double* BranchingBendingCosine::energy(double *coord, double *f, int *beadSet,
                                       double *kbend, double *eqt, double *z, int *params) {
        if(blocksnthreadse[1]>0) {
        BranchingBendingCosineenergy<<<blocksnthreadse[0], blocksnthreadse[1], (12 * blocksnthreadse[1]) * sizeof
                (double), stream>>> (coord, f, beadSet, kbend, eqt, params, gU_i, z, CUDAcommon::getCUDAvars()
                .gculpritID,
                CUDAcommon::getCUDAvars().gculpritFF,
                CUDAcommon::getCUDAvars().gculpritinteraction, gFF, ginteraction);
            CUDAcommon::handleerror( cudaGetLastError() ,"BranchingBendingCosineenergy", "BranchingBendingCosine.cu");
//        auto cvars = CUDAcommon::getCUDAvars();
//        cvars.streamvec.push_back(&stream);
//        CUDAcommon::cudavars = cvars;
//        double* gpu_Utot = CUDAcommon::getCUDAvars().gpu_energy;
//        addvector<<<1,1,0,stream>>>(gU_i,params, gU_sum, gpu_Utot);
//        CUDAcommon::handleerror( cudaGetLastError(),"BranchingBendingCosineenergy", "BranchingBendingCosine.cu");
//        return gU_sum;
        }
    if(blocksnthreadsez[1]>0) {

        BranchingBendingCosineenergyz << < blocksnthreadsez[0], blocksnthreadsez[1], (24 * blocksnthreadsez[1]) *
                                      sizeof(double), stream>> > (coord, f, beadSet, kbend, eqt, params, gU_i, z,
                CUDAcommon::getCUDAvars().gculpritID,
                CUDAcommon::getCUDAvars().gculpritFF,
                CUDAcommon::getCUDAvars().gculpritinteraction, gFF, ginteraction );
        CUDAcommon::handleerror(cudaGetLastError(),"BranchingBendingCosineenergyz", "BranchingBendingCosine.cu");
//        auto cvars = CUDAcommon::getCUDAvars();
//        cvars.streamvec.push_back(&stream);
//        CUDAcommon::cudavars = cvars;
//        double* gpu_Utot = CUDAcommon::getCUDAvars().gpu_energy;
//        addvector<<<1,1,0,stream>>>(gU_i,params, gU_sum, gpu_Utot);
//        CUDAcommon::handleerror(cudaGetLastError(),"BranchingBendingCosineenergyz", "BranchingBendingCosine.cu");
//        nvtxRangePushA("CCEBBzv2");
//        BranchingBendingCosineenergyz2 << < blocksnthreadsf[0], blocksnthreadsf[1], (12 * blocksnthreadsf[1]) *
//                                                                                     sizeof(double), stream>> > (coord, f, beadSet, kbend, eqt, params, gU_i, z );
//        CUDAcommon::handleerror(cudaGetLastError(),"BranchingBendingCosineenergyz2", "BranchingBendingCosine.cu");
//        nvtxRangePop();
//        return gU_sum;
    }
    if(blocksnthreadse[1]<=0 && blocksnthreadsez[1]<=0)
        return NULL;
    else {
        auto cvars = CUDAcommon::getCUDAvars();
        cvars.streamvec.push_back(&stream);
        CUDAcommon::cudavars = cvars;
        double *gpu_Utot = CUDAcommon::getCUDAvars().gpu_energy;

//        addvector<<<1,1,0,stream>>>(gU_i,params, gU_sum, gpu_Utot);
//        cudaStreamSynchronize(stream);
//        addvectorred<<<1,200,200*sizeof(double),stream>>>(gU_i,params, gU_sum, gpu_Utot);
//        cudaStreamSynchronize(stream);
//        std::cout<<"bntaddvec "<<bntaddvec2.at(0)<<" "<<bntaddvec2.at(1)<<" "<<bntaddvec2.at(0)<<" "
//                ""<<bntaddvec2.at(2)<<" "<<bntaddvec2.at(3)<<endl;
        resetdoublevariableCUDA << < 1, 1, 0, stream >> > (gU_sum);
        addvectorred2 << < bntaddvec2.at(2), bntaddvec2.at(3), bntaddvec2.at(3) * sizeof(double), stream >> > (gU_i,
                params, gU_sum, gpu_Utot);
//        CUDAcommon::handleerror(cudaDeviceSynchronize(),"FilamentBendingCosineenergyz", "FilamentBendingCosine.cu");
        CUDAcommon::handleerror(cudaGetLastError(), "FilamentBendingCosineenergyz", "FilamentBendingCosine.cu");
        return gU_sum;
    }
}

void BranchingBendingCosine::forces(double *coord, double *f, int *beadSet,
                                    double *kbend, double *eqt, int *params){
    if(blocksnthreadsf[1]>0) {
        BranchingBendingCosineforces << < blocksnthreadsf[0], blocksnthreadsf[1], (12 * blocksnthreadsf[1]) *
                                                                                  sizeof(double), stream >> > (coord, f, beadSet, kbend, eqt, params);
        auto cvars = CUDAcommon::getCUDAvars();
        cvars.streamvec.push_back(&stream);
        CUDAcommon::cudavars = cvars;
        CUDAcommon::handleerror(cudaGetLastError(),"BranchingBendingCosineforces", "BranchingBendingCosine.cu");
    }
}
void BranchingBendingCosine::checkforculprit() {
    CUDAcommon::printculprit("BranchingBending","BranchingBendingCosine");
    BranchingPoint* br;
    br = (BranchingPoint::getBranchingPoints()[CUDAcommon::getCUDAvars().culpritID[0]]);
    cout<<"Printing culprit branching point information."<<endl;
    br->printSelf();
    exit(EXIT_FAILURE);
}
#endif

double BranchingBendingCosine::energy(double *coord, double *f, int *beadSet,
                                      double *kbend, double *eqt){

    int n = BranchingBending<BranchingBendingCosine>::n;
    int nint = BranchingPoint::getBranchingPoints().size();

    double *coord1, *coord2, *coord3, *coord4, U_i, L1, L2, L1L2, l1l2, phi, dPhi;

    double U = 0.0;

    for(int i = 0; i < nint; i += 1) {

        coord1 = &coord[3 * beadSet[n * i]];
        coord2 = &coord[3 * beadSet[n * i + 1]];
        coord3 = &coord[3 * beadSet[n * i + 2]];
        coord4 = &coord[3 * beadSet[n * i + 3]];

        L1 = sqrt(scalarProduct(coord1, coord2,
                                coord1, coord2));
        L2 = sqrt(scalarProduct(coord3, coord4,
                                coord3, coord4));

        L1L2 = L1*L2;
        l1l2 = scalarProduct(coord1, coord2,
                             coord3, coord4);

        phi = safeacos(l1l2 / L1L2);
        dPhi = phi-eqt[i];

        U_i = kbend[i] * ( 1 - cos(dPhi) );

        if(fabs(U_i) == numeric_limits<double>::infinity()
           || U_i != U_i || U_i < -1.0) {

            //set culprit and return
            BranchingInteractions::_branchingCulprit = BranchingPoint::getBranchingPoints()[i];

            return -1;
        }

        U += U_i;
    }

    return U;
}

double BranchingBendingCosine::energy(double *coord, double *f, int *beadSet,
                                      double *kbend, double *eqt, double d){

    int n = BranchingBending<BranchingBendingCosine>::n;
    int nint = BranchingPoint::getBranchingPoints().size();

    double *coord1, *coord2, *coord3, *coord4, *force1, *force2, *force3, *force4, U_i, L1, L2, L1L2, l1l2, phi, dPhi;

    double U = 0.0;

    for(int i = 0; i < nint; i += 1) {

        coord1 = &coord[3 * beadSet[n * i]];
        coord2 = &coord[3 * beadSet[n * i + 1]];
        coord3 = &coord[3 * beadSet[n * i + 2]];
        coord4 = &coord[3 * beadSet[n * i + 3]];

        force1 = &f[3 * beadSet[n * i]];
        force2 = &f[3 * beadSet[n * i + 1]];
        force3 = &f[3 * beadSet[n * i + 2]];
        force4 = &f[3 * beadSet[n * i + 3]];


        L1 = sqrt(scalarProductStretched(coord1, force1, coord2, force2,
                                         coord1, force1, coord2, force2, d));
        L2 = sqrt(scalarProductStretched(coord3, force3, coord4, force4,
                                         coord3, force3, coord4, force4, d));

        L1L2 = L1*L2;
        l1l2 = scalarProductStretched(coord1, force1, coord2, force2,
                                      coord3, force3, coord4, force4, d);

        phi = safeacos(l1l2 / L1L2);
        dPhi = phi-eqt[i];

        U_i = kbend[i] * ( 1 - cos(dPhi) );

        if(fabs(U_i) == numeric_limits<double>::infinity()
           || U_i != U_i || U_i < -1.0) {

            //set culprit and return
            BranchingInteractions::_branchingCulprit = BranchingPoint::getBranchingPoints()[i];

            return -1;
        }

        U += U_i;
    }

    return U;
}

void BranchingBendingCosine::forces(double *coord, double *f, int *beadSet,
                                    double *kbend, double *eqt){


    int n = BranchingBending<BranchingBendingCosine>::n;
    int nint = BranchingPoint::getBranchingPoints().size();

    double *coord1, *coord2, *coord3, *coord4, *force1, *force2, *force3, *force4;
    double L1, L2, L1L2, l1l2, phi, dPhi, A, B, C, invL1, invL2, k;

//    double U = 0;

    for(int i = 0; i < nint; i += 1) {

        coord1 = &coord[3 * beadSet[n * i]];
        coord2 = &coord[3 * beadSet[n * i + 1]];
        coord3 = &coord[3 * beadSet[n * i + 2]];
        coord4 = &coord[3 * beadSet[n * i + 3]];

        force1 = &f[3 * beadSet[n * i]];
        force2 = &f[3 * beadSet[n * i + 1]];
        force3 = &f[3 * beadSet[n * i + 2]];
        force4 = &f[3 * beadSet[n * i + 3]];

        L1 = sqrt(scalarProduct(coord1, coord2,
                                coord1, coord2));
        L2 = sqrt(scalarProduct(coord3, coord4,
                                coord3, coord4));

        L1L2 = L1*L2;
        l1l2 = scalarProduct(coord1, coord2,
                             coord3, coord4);

        invL1 = 1/L1;
        invL2 = 1/L2;
        A = invL1*invL2;
        B = l1l2*invL1*A*A*L2;
        C = l1l2*invL2*A*A*L1;

//        phi = safeacos(l1l2 / L1L2);
        phi = safeacos(l1l2 * A);
        dPhi = phi-eqt[i];

        k =  kbend[i] * sin(dPhi)/sin(phi);

        //force on i, f = k*(-A*l2 + 2*B*l1):
        force1[0] += k * ((coord3[0] - coord4[0])*A +
                          (coord2[0] - coord1[0])*B );
        force1[1] += k * ((coord3[1] - coord4[1])*A +
                          (coord2[1] - coord1[1])*B );
        force1[2] += k * ((coord3[2] - coord4[2])*A +
                          (coord2[2] - coord1[2])*B );


        //force on i+1, f = k*(A*l2 - 2*B*l1):
        force2[0] += k * ((-coord3[0] + coord4[0])*A -
                          (coord2[0] - coord1[0])*B );
        force2[1] += k * ((-coord3[1] + coord4[1])*A -
                          (coord2[1] - coord1[1])*B );
        force2[2] += k * ((-coord3[2] + coord4[2])*A -
                          (coord2[2] - coord1[2])*B );

        //force on j, k*(-A*l1 + 2*C*l2):
        force3[0] += k *((coord1[0] - coord2[0])*A +
                         (coord4[0] - coord3[0])*C );
        force3[1] += k *((coord1[1] - coord2[1])*A +
                         (coord4[1] - coord3[1])*C );
        force3[2] += k *((coord1[2] - coord2[2])*A +
                         (coord4[2] - coord3[2])*C );

        //force on j+1, k*(A*l1 - 2*C*l2):
        force4[0] += k *((-coord1[0] + coord2[0])*A -
                         (coord4[0] - coord3[0])*C );
        force4[1] += k *((-coord1[1] + coord2[1])*A -
                         (coord4[1] - coord3[1])*C );
        force4[2] += k *((-coord1[2] + coord2[2])*A -
                         (coord4[2] - coord3[2])*C );
    }
}
