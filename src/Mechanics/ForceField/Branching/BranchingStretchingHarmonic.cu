
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

#include "BranchingStretchingHarmonic.h"
#include "BranchingStretchingHarmonicCUDA.h"
#include "BranchingStretching.h"

#include "BranchingPoint.h"
#include "Bead.h"

#include "MathFunctions.h"
#include <cuda.h>
#include <cuda_runtime.h>
#include "nvToolsExt.h"

using namespace mathfunc;
#ifdef CUDAACCL
void BranchingStretchingHarmonic::deallocate(){
    if(!(CUDAcommon::getCUDAvars().conservestreams))
        CUDAcommon::handleerror(cudaStreamDestroy(stream));
    CUDAcommon::handleerror(cudaFree(gU_i));
    CUDAcommon::handleerror(cudaFree(gU_sum));
    CUDAcommon::handleerror(cudaFree(gFF));
    CUDAcommon::handleerror(cudaFree(ginteraction));
}
void BranchingStretchingHarmonic::optimalblocksnthreads( int nint){
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
                                                       BranchingStretchingHarmonicenergy, blockToSmemFB, 0);
        blocksnthreadse.push_back((nint + blockSize - 1) / blockSize);
        blocksnthreadse.push_back(blockSize);
        blockSize = 0;

        cudaOccupancyMaxPotentialBlockSizeVariableSMem(&minGridSize, &blockSize,
                                                       BranchingStretchingHarmonicenergyz, blockToSmemFB2, 0);
        blocksnthreadsez.push_back((nint + blockSize - 1) / blockSize);
        blocksnthreadsez.push_back(blockSize);
        blockSize = 0;

        cudaOccupancyMaxPotentialBlockSizeVariableSMem(&minGridSize, &blockSize,
                                                       BranchingStretchingHarmonicforces, blockToSmemFB, 0);
        blocksnthreadsf.push_back((nint + blockSize - 1) / blockSize);
        blocksnthreadsf.push_back(blockSize);
        //get addition vars
        bntaddvec2.clear();
        bntaddvec2 = getaddred2bnt(nint);
        CUDAcommon::handleerror(cudaMalloc((void **) &gU_i, bntaddvec2.at(0)*sizeof(double)));
        CUDAcommon::handleerror(cudaMemset(gU_i, 0, bntaddvec2.at(0) * sizeof(double)));
//        CUDAcommon::handleerror(cudaMalloc((void **) &gU_i, nint*sizeof(double)));
        CUDAcommon::handleerror(cudaMalloc((void **) &gU_sum, sizeof(double)));

//        CUDAcommon::handleerror(cudaMalloc((void **) &gU_i, nint*sizeof(double)));
//        CUDAcommon::handleerror(cudaMalloc((void **) &gU_sum, sizeof(double)));

        char a[] = "BranchingFF";
        char b[] = "Branching Stretching Harmonic";
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
double* BranchingStretchingHarmonic::energy(double *coord, double *f, int *beadSet,
                                            double *kstr, double *eql, double *pos, int *params) {
//    if(blocksnthreadse[1]>0) {
//        BranchingStretchingHarmonicenergy<<<blocksnthreadse[0], blocksnthreadse[1], (9 * blocksnthreadse[1]) * sizeof
//                (double), stream>>> (coord, f, beadSet, kstr, eql, pos, params, gU_i, CUDAcommon::getCUDAvars()
//                .gculpritID, CUDAcommon::getCUDAvars().gculpritFF, CUDAcommon::getCUDAvars().gculpritinteraction,
//                gFF, ginteraction);
//        auto cvars = CUDAcommon::getCUDAvars();
//        cvars.streamvec.push_back(&stream);
//        CUDAcommon::cudavars = cvars;
//        CUDAcommon::handleerror( cudaGetLastError() ,"BranchingStretchingHarmonicenergy", "BranchingStretchingHarmonic.cu");
//        double* gpu_Utot = CUDAcommon::getCUDAvars().gpu_energy;
//        addvector<<<1,1,0,stream>>>(gU_i,params, gU_sum, gpu_Utot);
//        CUDAcommon::handleerror( cudaGetLastError() ,"BranchingStretchingHarmonicenergy", "BranchingStretchingHarmonic.cu");
//        return gU_sum;}
//    else
//        return NULL;
}


double* BranchingStretchingHarmonic::energy(double *coord, double *f, int *beadSet,
                                            double *kstr, double *eql, double *pos, double *z, int *params) {
    if(blocksnthreadse[1]>0) {
        BranchingStretchingHarmonicenergy<<<blocksnthreadse[0], blocksnthreadse[1], (9 * blocksnthreadse[1]) * sizeof
                (double), stream>>> (coord, f, beadSet, kstr, eql, pos, params, gU_i, z, CUDAcommon::getCUDAvars()
                .gculpritID, CUDAcommon::getCUDAvars().gculpritFF, CUDAcommon::getCUDAvars().gculpritinteraction,
                gFF, ginteraction);
//        auto cvars = CUDAcommon::getCUDAvars();
//        cvars.streamvec.push_back(&stream);
//        CUDAcommon::cudavars = cvars;
//        CUDAcommon::handleerror( cudaGetLastError() ,"BranchingStretchingHarmonicenergy", "BranchingStretchingHarmonic.cu");
//        double* gpu_Utot = CUDAcommon::getCUDAvars().gpu_energy;
//        addvector<<<1,1,0,stream>>>(gU_i,params, gU_sum, gpu_Utot);
//        CUDAcommon::handleerror( cudaGetLastError() ,"BranchingStretchingHarmonicenergy", "BranchingStretchingHarmonic.cu");
//        return gU_sum;
    }

    if(blocksnthreadsez[1]>0) {
        BranchingStretchingHarmonicenergyz << < blocksnthreadsez[0], blocksnthreadsez[1], (18 * blocksnthreadsez[1]) *
                                                                                          sizeof(double), stream>> >
                (coord, f, beadSet, kstr, eql, pos, params, gU_i, z, CUDAcommon::getCUDAvars().gculpritID,
                CUDAcommon::getCUDAvars().gculpritFF,
                CUDAcommon::getCUDAvars().gculpritinteraction, gFF, ginteraction);
//        auto cvars = CUDAcommon::getCUDAvars();
//        cvars.streamvec.push_back(&stream);
//        CUDAcommon::cudavars = cvars;
//        CUDAcommon::handleerror(cudaGetLastError(),"BranchingStretchingHarmonicenergyz", "BranchingStretchingHarmonic"
//                ".cu");
//        double* gpu_Utot = CUDAcommon::getCUDAvars().gpu_energy;
//
//        addvector<<<1,1,0,stream>>>(gU_i,params, gU_sum, gpu_Utot);
//        CUDAcommon::handleerror(cudaGetLastError(),"BranchingStretchingHarmonicenergyz", "BranchingStretchingHarmonic"
//                ".cu");
//        return gU_sum;
    }
    if(blocksnthreadse[1]<=0 && blocksnthreadsez[1]<=0)
        return NULL;
    else{
        auto cvars = CUDAcommon::getCUDAvars();
        cvars.streamvec.push_back(&stream);
        CUDAcommon::cudavars = cvars;
        double* gpu_Utot = CUDAcommon::getCUDAvars().gpu_energy;

//        addvector<<<1,1,0,stream>>>(gU_i,params, gU_sum, gpu_Utot);
//        cudaStreamSynchronize(stream);
//        addvectorred<<<1,200,200*sizeof(double),stream>>>(gU_i,params, gU_sum, gpu_Utot);
//        cudaStreamSynchronize(stream);
//        std::cout<<"bntaddvec "<<bntaddvec2.at(0)<<" "<<bntaddvec2.at(1)<<" "<<bntaddvec2.at(0)<<" "
//                ""<<bntaddvec2.at(2)<<" "<<bntaddvec2.at(3)<<endl;
        resetdoublevariableCUDA<<<1,1,0,stream>>>(gU_sum);
        addvectorred2<<<bntaddvec2.at(2),bntaddvec2.at(3), bntaddvec2.at(3) * sizeof(double),stream>>>(gU_i,
                params, gU_sum, gpu_Utot);
//        CUDAcommon::handleerror(cudaDeviceSynchronize(),"FilamentBendingCosineenergyz", "FilamentBendingCosine.cu");
        CUDAcommon::handleerror(cudaGetLastError(),"FilamentBendingCosineenergyz", "FilamentBendingCosine.cu");
        return gU_sum;
    }

}

void BranchingStretchingHarmonic::forces(double *coord, double *f, int *beadSet,
                                         double *kstr, double *eql, double *pos, int *params){
    if(blocksnthreadsf[1]>0) {
        BranchingStretchingHarmonicforces << < blocksnthreadsf[0], blocksnthreadsf[1], (9 * blocksnthreadsf[1]) *
                                                                                       sizeof(double), stream >> >
                (coord, f, beadSet, kstr, eql, pos, params);
        auto cvars = CUDAcommon::getCUDAvars();
        cvars.streamvec.push_back(&stream);
        CUDAcommon::cudavars = cvars;
    }
}
void BranchingStretchingHarmonic::checkforculprit() {
    CUDAcommon::printculprit("BranchingStretching","BranchingStretchingHarmonic");
    BranchingPoint* br;
    br = (BranchingPoint::getBranchingPoints()[CUDAcommon::getCUDAvars().culpritID[0]]);
    cout<<"Printing culprit branching point information."<<endl;
    br->printSelf();
    exit(EXIT_FAILURE);
}
#endif
double BranchingStretchingHarmonic::energy(double *coord, double *f, int *beadSet,
                                           double *kstr, double *eql, double *pos){

    int n = BranchingStretching<BranchingStretchingHarmonic>::n;
    int nint = BranchingPoint::getBranchingPoints().size();

    double *coord1, *coord2, *coord3, dist, U_i;
    double *v1 = new double[3];

    double U = 0;

    for(int i = 0; i < nint; i += 1) {

        coord1 = &coord[3 * beadSet[n * i]];
        coord2 = &coord[3 * beadSet[n * i + 1]];
        coord3 = &coord[3 * beadSet[n * i + 2]];

        midPointCoordinate(v1, coord1, coord2, pos[i]);
        dist = twoPointDistance(v1, coord3) - eql[i];

        U_i = 0.5 * kstr[i] * dist * dist;


        if(fabs(U_i) == numeric_limits<double>::infinity()
           || U_i != U_i || U_i < -1.0) {

            //set culprit and return
            BranchingInteractions::_branchingCulprit = BranchingPoint::getBranchingPoints()[i];

            return -1;
        }

        U += U_i;
    }
    delete v1;
    return U;
}

double BranchingStretchingHarmonic::energy(double *coord, double *f, int *beadSet,
                                           double *kstr, double *eql, double *pos, double d){

    int n = BranchingStretching<BranchingStretchingHarmonic>::n;
    int nint = BranchingPoint::getBranchingPoints().size();

    double *coord1, *coord2, *coord3, *f1, *f2, *f3, dist, U_i;
    double *v1 = new double[3];
    double *vzero = new double[3]; vzero[0] = 0; vzero[1] = 0; vzero[2] = 0;

    double U = 0;

    for(int i = 0; i < nint; i += 1) {

        coord1 = &coord[3 * beadSet[n * i]];
        coord2 = &coord[3 * beadSet[n * i + 1]];
        coord3 = &coord[3 * beadSet[n * i + 2]];

        f1 = &f[3 * beadSet[n * i]];
        f2 = &f[3 * beadSet[n * i + 1]];
        f3 = &f[3 * beadSet[n * i + 2]];

        midPointCoordinateStretched(v1, coord1, f1, coord2, f2, pos[i], d);
        dist = twoPointDistanceStretched(v1, vzero, coord3, f3, d) - eql[i];

        U_i = 0.5 * kstr[i] * dist * dist;


        if(fabs(U_i) == numeric_limits<double>::infinity()
           || U_i != U_i || U_i < -1.0) {

            //set culprit and return
            BranchingInteractions::_branchingCulprit = BranchingPoint::getBranchingPoints()[i];

            return -1;
        }

        U += U_i;
    }
    delete v1; vzero;
    return U;
}

void BranchingStretchingHarmonic::forces(double *coord, double *f, int *beadSet,
                                         double *kstr, double *eql, double *pos,
                                         double *stretchforce){


    int n = BranchingStretching<BranchingStretchingHarmonic>::n;
    int nint = BranchingPoint::getBranchingPoints().size();

    double *coord1, *coord2, *coord3, *f1, *f2, *f3, dist, invL, f0;
    double *v1 = new double[3];


    for(int i = 0; i < nint; i += 1) {

        coord1 = &coord[3 * beadSet[n * i]];
        coord2 = &coord[3 * beadSet[n * i + 1]];
        coord3 = &coord[3 * beadSet[n * i + 2]];

        f1 = &f[3 * beadSet[n * i]];
        f2 = &f[3 * beadSet[n * i + 1]];
        f3 = &f[3 * beadSet[n * i + 2]];


        midPointCoordinate(v1, coord1, coord2, pos[i]);
        dist = twoPointDistance(v1, coord3);

        invL = 1 / dist;
        f0 = kstr[i] * ( dist - eql[i]) * invL;

        f1[0] +=  -f0 * ( coord3[0] - v1[0] ) * (pos[i] - 1);
        f1[1] +=  -f0 * ( coord3[1] - v1[1] ) * (pos[i] - 1);
        f1[2] +=  -f0 * ( coord3[2] - v1[2] ) * (pos[i] - 1);

        // force i+1
        f2[0] +=  f0 * ( coord3[0] - v1[0] ) * pos[i];
        f2[1] +=  f0 * ( coord3[1] - v1[1] ) * pos[i];
        f2[2] +=  f0 * ( coord3[2] - v1[2] ) * pos[i];

        //force on j
        f3[0] +=  -f0 * ( coord3[0] - v1[0] );
        f3[1] +=  -f0 * ( coord3[1] - v1[1] );
        f3[2] +=  -f0 * ( coord3[2] - v1[2] );
        stretchforce[i] = f0/invL;

    }
    delete v1;
}
