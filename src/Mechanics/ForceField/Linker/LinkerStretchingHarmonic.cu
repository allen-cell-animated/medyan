
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

#include "LinkerStretchingHarmonic.h"
#include "LinkerStretchingHarmonicCUDA.h"
#include "LinkerStretching.h"
#include "Linker.h"

#include "Bead.h"

#include "MathFunctions.h"
#include <cuda.h>
#include <cuda_runtime.h>
#include "nvToolsExt.h"

using namespace mathfunc;
#ifdef CUDAACCL
void LinkerStretchingHarmonic::deallocate(){
    CUDAcommon::handleerror(cudaStreamDestroy(stream));
    CUDAcommon::handleerror(cudaFree(gU_i));
    CUDAcommon::handleerror(cudaFree(gU_sum));
    CUDAcommon::handleerror(cudaFree(gFF));
    CUDAcommon::handleerror(cudaFree(ginteraction));
}
void LinkerStretchingHarmonic::optimalblocksnthreads( int nint){
    //CUDA stream create
    CUDAcommon::handleerror(cudaStreamCreate(&stream));
    blocksnthreadse.clear();
    blocksnthreadsez.clear();
    blocksnthreadsf.clear();
    int blockSize;   // The launch configurator returned block size
    int minGridSize; // The minimum grid size needed to achieve the
    // maximum occupancy for a full device launch
//    int gridSize;    // The actual grid size needed, based on input size
//    unaryfn::argument_type blksize;
//    unaryfn::result_type result;
//    unaryfn ufn;
    if(nint>0) {
        numint = nint;
        cudaOccupancyMaxPotentialBlockSizeVariableSMem(&minGridSize, &blockSize,
                                                       LinkerStretchingHarmonicenergy, blockToSmem, 0);
        blocksnthreadse.push_back((nint + blockSize - 1) / blockSize);
        blocksnthreadse.push_back(blockSize);
        blockSize = 0;

        cudaOccupancyMaxPotentialBlockSizeVariableSMem(&minGridSize, &blockSize,
                                                       LinkerStretchingHarmonicenergyz, blockToSmemez, 0);
        blocksnthreadsez.push_back((nint + blockSize - 1) / blockSize);
        blocksnthreadsez.push_back(blockSize);
        blockSize = 0;

        cudaOccupancyMaxPotentialBlockSizeVariableSMem(&minGridSize, &blockSize,
                                                       LinkerStretchingHarmonicforces, blockToSmem, 0);
        blocksnthreadsf.push_back((nint + blockSize - 1) / blockSize);
        blocksnthreadsf.push_back(blockSize);
        //get addition vars
        bntaddvec2.clear();
        bntaddvec2 = getaddred2bnt(nint);
        CUDAcommon::handleerror(cudaMalloc((void **) &gU_i, bntaddvec2.at(0)*sizeof(double)));
        CUDAcommon::handleerror(cudaMemset(gU_i, 0, bntaddvec2.at(0) * sizeof(double)));
//        CUDAcommon::handleerror(cudaMalloc((void **) &gU_i, nint*sizeof(double)));
        CUDAcommon::handleerror(cudaMalloc((void **) &gU_sum, sizeof(double)));
        char a[] = "LinkerFF";
        char b[] = "Linker Stretching Harmonic";
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
double* LinkerStretchingHarmonic::energy(double *coord, double *f, int *beadSet,
                                         double *kstr, double *eql, double *pos1, double *pos2,
                                         int *params) {
//    if(blocksnthreadse[1]>0) {
//
//        LinkerStretchingHarmonicenergy<<<blocksnthreadse[0], blocksnthreadse[1], (12 * blocksnthreadse[1]) * sizeof
//                (double), stream>>>
//                          (coord, f, beadSet, kstr, eql, pos1, pos2, params, gU_i, CUDAcommon::getCUDAvars().gculpritID,
//                                  CUDAcommon::getCUDAvars().gculpritFF,
//                                  CUDAcommon::getCUDAvars().gculpritinteraction, gFF, ginteraction);
//        CUDAcommon::handleerror( cudaGetLastError() ,"LinkerStretchingHarmonicenergy", "LinkerStretchingHarmonic.cu");
//        auto cvars = CUDAcommon::getCUDAvars();
//        cvars.streamvec.push_back(&stream);
//        CUDAcommon::cudavars = cvars;
//        double* gpu_Utot = CUDAcommon::getCUDAvars().gpu_energy;
//        addvector<<<1,1,0,stream>>>(gU_i,params, gU_sum, gpu_Utot);
//        CUDAcommon::handleerror( cudaGetLastError() ,"LinkerStretchingHarmonicenergy", "LinkerStretchingHarmonic.cu");
//        return gU_sum;}
//    else
//        return NULL;
}


double* LinkerStretchingHarmonic::energy(double *coord, double *f, int *beadSet,
                                         double *kstr, double *eql, double *pos1, double *pos2, double *z,
                                         int *params) {
//    nvtxRangePushA("E_wait");
//    CUDAcommon::handleerror(cudaStreamWaitEvent(stream, *(CUDAcommon::getCUDAvars().event), 0));
//    nvtxRangePop();
    if(blocksnthreadse[1]>0) {

        LinkerStretchingHarmonicenergy<<<blocksnthreadse[0], blocksnthreadse[1], (12 * blocksnthreadse[1]) * sizeof
                (double), stream>>>
                          (coord, f, beadSet, kstr, eql, pos1, pos2, params, gU_i, z, CUDAcommon::getCUDAvars()
                                  .gculpritID,
                                  CUDAcommon::getCUDAvars().gculpritFF,
                                  CUDAcommon::getCUDAvars().gculpritinteraction, gFF, ginteraction);
        CUDAcommon::handleerror( cudaGetLastError() ,"LinkerStretchingHarmonicenergy", "LinkerStretchingHarmonic.cu");
//        auto cvars = CUDAcommon::getCUDAvars();
//        cvars.streamvec.push_back(&stream);
//        CUDAcommon::cudavars = cvars;
//        double* gpu_Utot = CUDAcommon::getCUDAvars().gpu_energy;
//        addvector<<<1,1,0,stream>>>(gU_i,params, gU_sum, gpu_Utot);
//        CUDAcommon::handleerror( cudaGetLastError() ,"LinkerStretchingHarmonicenergy", "LinkerStretchingHarmonic.cu");
//        return gU_sum;
    }

    if(blocksnthreadsez[1]>0) {
        LinkerStretchingHarmonicenergyz << < blocksnthreadsez[0], blocksnthreadsez[1], (24 * blocksnthreadsez[1]) *
                                             sizeof(double), stream>> >
                                            (coord, f, beadSet, kstr, eql, pos1, pos2, params, gU_i, z ,
                                             CUDAcommon::getCUDAvars().gculpritID,
                                             CUDAcommon::getCUDAvars().gculpritFF,
                                             CUDAcommon::getCUDAvars().gculpritinteraction, gFF, ginteraction);
        CUDAcommon::handleerror(cudaGetLastError(),"LinkerStretchingHarmonicenergyz", "LinkerStretchingHarmonic.cu");
//        auto cvars = CUDAcommon::getCUDAvars();
//        cvars.streamvec.push_back(&stream);
//        CUDAcommon::cudavars = cvars;
//
//        double* gpu_Utot = CUDAcommon::getCUDAvars().gpu_energy;
//        addvector<<<1,1,0,stream>>>(gU_i,params, gU_sum, gpu_Utot);
//        CUDAcommon::handleerror(cudaGetLastError(),"LinkerStretchingHarmonicenergyz", "LinkerStretchingHarmonic.cu");
//
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
//        std::cout<<"bntaddvec "<<bntaddvec2.at(0)<<" "<<bntaddvec2.at(1)<<" "<<bntaddvec2.at(2)<<" "<<bntaddvec2.at
//                (3)<<" "<<numint<<endl;
        resetdoublevariableCUDA<<<1,1,0,stream>>>(gU_sum);
        addvectorred2<<<bntaddvec2.at(2),bntaddvec2.at(3), bntaddvec2.at(3) * sizeof(double),stream>>>(gU_i,
                params, gU_sum, gpu_Utot);
        CUDAcommon::handleerror( cudaGetLastError() ,"LinkerStretchingHarmonicenergy", "LinkerStretchingHarmonic.cu");
        return gU_sum;
    }
}
void LinkerStretchingHarmonic::forces(double *coord, double *f, int *beadSet,
                                      double *kstr, double *eql, double *pos1, double *pos2, int *params) {
    if (blocksnthreadsf[1] > 0) {
        LinkerStretchingHarmonicforces << < blocksnthreadsf[0], blocksnthreadsf[1], (12 *
                                                                                     blocksnthreadsf[1]) *
                                                                                    sizeof(double), stream >> >
                                                                                                    (coord, f, beadSet, kstr, eql, pos1, pos2, params);
        auto cvars = CUDAcommon::getCUDAvars();
        cvars.streamvec.push_back(&stream);
        CUDAcommon::cudavars = cvars;
        CUDAcommon::handleerror(cudaGetLastError(),"LinkerStretchingHarmonicforces", "LinkerStretchingHarmonic.cu");
    }
}
void LinkerStretchingHarmonic::checkforculprit() {
    CUDAcommon::printculprit("LinkerStretching","LinkerStretchingHarmonic");
    Linker* l;
    l = Linker::getLinkers()[CUDAcommon::getCUDAvars().culpritID[0]];
    cout<<"Printing culprit Filament information."<<endl;
    l->printSelf();
    exit(EXIT_FAILURE);
}

#endif
    double LinkerStretchingHarmonic::energy(double *coord, double *f, int *beadSet,
                                            double *kstr, double *eql, double *pos1, double *pos2) {

        int n = LinkerStretching<LinkerStretchingHarmonic>::n;
        int nint = Linker::getLinkers().size();

        double *coord1, *coord2, *coord3, *coord4, dist, U_i;
        double *v1 = new double[3];
        double *v2 = new double[3];

        double U = 0;

        for(int i = 0; i < nint; i += 1) {

            coord1 = &coord[3 * beadSet[n * i]];
            coord2 = &coord[3 * beadSet[n * i + 1]];
            coord3 = &coord[3 * beadSet[n * i + 2]];
            coord4 = &coord[3 * beadSet[n * i + 3]];

            midPointCoordinate(v1, coord1, coord2, pos1[i]);
            midPointCoordinate(v2, coord3, coord4, pos2[i]);

            dist = twoPointDistance(v1, v2) - eql[i];
            U_i = 0.5 * kstr[i] * dist * dist;

            if(fabs(U_i) == numeric_limits<double>::infinity()
               || U_i != U_i || U_i < -1.0) {

                //set culprit and return
                LinkerInteractions::_linkerCulprit = Linker::getLinkers()[i];

                return -1;
            }

            U += U_i;
        }
        delete v1;
        delete v2;

        return U;
    }

    double LinkerStretchingHarmonic::energy(double *coord, double * f, int *beadSet,
                                            double *kstr, double *eql, double *pos1, double *pos2, double d){

        int n = LinkerStretching<LinkerStretchingHarmonic>::n;
        int nint = Linker::getLinkers().size();

        double *coord1, *coord2, *coord3, *coord4, *f1, *f2, *f3, *f4, dist;
        double *v1 = new double[3];
        double *v2 = new double[3];

        double U = 0;

        for(int i = 0; i < nint; i += 1) {

            coord1 = &coord[3 * beadSet[n * i]];
            coord2 = &coord[3 * beadSet[n * i + 1]];
            coord3 = &coord[3 * beadSet[n * i + 2]];
            coord4 = &coord[3 * beadSet[n * i + 3]];

            f1 = &f[3 * beadSet[n * i]];
            f2 = &f[3 * beadSet[n * i + 1]];
            f3 = &f[3 * beadSet[n * i + 2]];
            f4 = &f[3 * beadSet[n * i + 3]];

            midPointCoordinateStretched(v1, coord1, f1, coord2, f2, pos1[i], d);
            midPointCoordinateStretched(v2, coord3, f3, coord4, f4, pos2[i], d);

            dist = twoPointDistance(v1,  v2) - eql[i];

            U += 0.5 * kstr[i] * dist * dist;
        }
        delete v1;
        delete v2;

        return U;

    }
    void LinkerStretchingHarmonic::forces(double *coord, double *f, int *beadSet,
                                          double *kstr, double *eql, double *pos1, double *pos2){


        int n = LinkerStretching<LinkerStretchingHarmonic>::n;
        int nint = Linker::getLinkers().size();

        double *coord1, *coord2, *coord3, *coord4, dist, invL;
        double *v1 = new double[3];
        double *v2 = new double[3];

        double f0, *f1, *f2, *f3, *f4;

        for(int i = 0; i < nint; i += 1) {

            coord1 = &coord[3 * beadSet[n * i]];
            coord2 = &coord[3 * beadSet[n * i + 1]];
            coord3 = &coord[3 * beadSet[n * i + 2]];
            coord4 = &coord[3 * beadSet[n * i + 3]];

            midPointCoordinate(v1, coord1, coord2, pos1[i]);
            midPointCoordinate(v2, coord3, coord4, pos2[i]);

            dist = twoPointDistance(v1, v2) ;
            invL = 1 / dist;

            f0 = kstr[i] * ( dist - eql[i] ) * invL;

            f1 = &f[3 * beadSet[n * i]];
            f2 = &f[3 * beadSet[n * i + 1]];
            f3 = &f[3 * beadSet[n * i + 2]];
            f4 = &f[3 * beadSet[n * i + 3]];

            //force on i
            f1[0] +=   -f0 * ( v1[0] - v2[0] ) * (1 - pos1[i]);
            f1[1] +=   -f0 * ( v1[1] - v2[1] ) * (1 - pos1[i]);
            f1[2] +=   -f0 * ( v1[2] - v2[2] ) * (1 - pos1[i]);

            // force i+1
            f2[0] +=   -f0 * ( v1[0] - v2[0] ) * (pos1[i]);
            f2[1] +=   -f0 * ( v1[1] - v2[1] ) * (pos1[i]);
            f2[2] +=   -f0 * ( v1[2] - v2[2] ) * (pos1[i]);

            //force on j
            f3[0] +=   f0 * ( v1[0] - v2[0] ) * (1 - pos2[i]);
            f3[1] +=   f0 * ( v1[1] - v2[1] ) * (1 - pos2[i]);
            f3[2] +=   f0 * ( v1[2] - v2[2] ) * (1 - pos2[i]);

            // force j+1
            f4[0] +=   f0 * ( v1[0] - v2[0] ) * (pos2[i]);
            f4[1] +=   f0 * ( v1[1] - v2[1] ) * (pos2[i]);
            f4[2] +=   f0 * ( v1[2] - v2[2] ) * (pos2[i]);

            //assign stretch force
            Linker::getLinkers()[i]->getMLinker()->stretchForce = f0;
//        std::cout<<"LINKER "<<f1[0]<<" "<<f1[1]<<" "<<f1[2]<<" "<<f2[0]<<" "<<f2[1]<<" "<<f2[2]<<" "<<f3[0]<<" "<<f3[1]<<" "<<f3[2]<<" "<<f4[0]<<" "<<f4[1]<<" "<<f4[2]<<endl;
        }
        delete v1;
        delete v2;
    }

