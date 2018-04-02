
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

#include "FilamentStretchingHarmonic.h"
#include "FilamentStretching.h"
#include "FilamentStretchingHarmonicCUDA.h"
#include "Cylinder.h"
#include "Filament.h"
#include "Bead.h"

#include "MathFunctions.h"
#include <cuda.h>
#include <cuda_runtime.h>
#include "nvToolsExt.h"

using namespace mathfunc;
#ifdef CUDAACCL
void FilamentStretchingHarmonic::deallocate(){
    CUDAcommon::handleerror(cudaStreamDestroy(stream));
    CUDAcommon::handleerror(cudaFree(gU_i));
    CUDAcommon::handleerror(cudaFree(gU_sum));
    CUDAcommon::handleerror(cudaFree(gFF));
    CUDAcommon::handleerror(cudaFree(ginteraction));
}
void FilamentStretchingHarmonic::optimalblocksnthreads( int nint){
    //CUDA stream create
    CUDAcommon::handleerror(cudaStreamCreate(&stream));
    blocksnthreadse.clear();
    blocksnthreadsez.clear();
    blocksnthreadsf.clear();
    int blockSize;   // The launch configurator returned block size
    int minGridSize; // The minimum grid size needed to achieve the
    // maximum occupancy for a full device launch
    if(nint>0) {
        cudaOccupancyMaxPotentialBlockSizeVariableSMem(&minGridSize, &blockSize,
                                                       FilamentStretchingHarmonicenergy, blockToSmemF, 0);
        blocksnthreadse.push_back((nint + blockSize - 1) / blockSize);
        blocksnthreadse.push_back(blockSize);
        blockSize = 0;

        cudaOccupancyMaxPotentialBlockSizeVariableSMem(&minGridSize, &blockSize,
                                                       FilamentStretchingHarmonicenergyz, blockToSmem, 0);
        blocksnthreadsez.push_back((nint + blockSize - 1) / blockSize);
        blocksnthreadsez.push_back(blockSize);
        blockSize = 0;

        cudaOccupancyMaxPotentialBlockSizeVariableSMem(&minGridSize, &blockSize,
                                                       FilamentStretchingHarmonicforces, blockToSmemF, 0);
        blocksnthreadsf.push_back((nint + blockSize - 1) / blockSize);
        blocksnthreadsf.push_back(blockSize);
        //get addition vars
        bntaddvec2.clear();
        bntaddvec2 = getaddred2bnt(nint);
        CUDAcommon::handleerror(cudaMalloc((void **) &gU_i, bntaddvec2.at(0)*sizeof(double)));
        CUDAcommon::handleerror(cudaMemset(gU_i, 0, bntaddvec2.at(0) * sizeof(double)));
        CUDAcommon::handleerror(cudaMalloc((void **) &gU_sum, sizeof(double)));
        char a[] = "FilamentFF";
        char b[] = "Filament Stretching Harmonic";
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
double* FilamentStretchingHarmonic::energy(double *coord, double *f, int *beadSet,
                                             double *kstr, double *eql, int *params) {
//    if(blocksnthreadse[1]>0) {
//        FilamentStretchingHarmonicenergy<<<blocksnthreadse[0], blocksnthreadse[1], (6 * blocksnthreadse[1]) * sizeof
//                (double), stream>>> (coord, f, beadSet, kstr, eql, params, gU_i, CUDAcommon::getCUDAvars()
//                .gculpritID,
//                CUDAcommon::getCUDAvars().gculpritFF,
//                CUDAcommon::getCUDAvars().gculpritinteraction, gFF, ginteraction);
//        auto cvars = CUDAcommon::getCUDAvars();
//        cvars.streamvec.push_back(&stream);
//        CUDAcommon::cudavars = cvars;
//        CUDAcommon::handleerror( cudaGetLastError(),"FilamentStretchingHarmonicenergy", "FilamentStretchingHarmonic.cu");
//        double* gpu_Utot = CUDAcommon::getCUDAvars().gpu_energy;
//        addvector<<<1,1,0,stream>>>(gU_i, params, gU_sum, gpu_Utot);
//        CUDAcommon::handleerror( cudaGetLastError() ,"FilamentStretchingHarmonicenergy", "FilamentStretchingHarmonic.cu");
//
//        return gU_sum;
//    }
//    else
//        return NULL;
}


double* FilamentStretchingHarmonic::energy(double *coord, double *f, int *beadSet,
                                             double *kstr, double *eql, double *z, int *params) {
//    nvtxRangePushA("E_wait");
//    CUDAcommon::handleerror(cudaStreamWaitEvent(stream, *(CUDAcommon::getCUDAvars().event), 0));
//    nvtxRangePop();
    if(blocksnthreadse[1]>0) {
        FilamentStretchingHarmonicenergy<<<blocksnthreadse[0], blocksnthreadse[1], (6 * blocksnthreadse[1]) * sizeof
                (double), stream>>> (coord, f, beadSet, kstr, eql, params, gU_i, z, CUDAcommon::getCUDAvars()
                .gculpritID,
                CUDAcommon::getCUDAvars().gculpritFF,
                CUDAcommon::getCUDAvars().gculpritinteraction, gFF, ginteraction);
        CUDAcommon::handleerror( cudaGetLastError(),"FilamentStretchingHarmonicenergy", "FilamentStretchingHarmonic.cu");
//        auto cvars = CUDAcommon::getCUDAvars();
//        cvars.streamvec.push_back(&stream);
//        CUDAcommon::cudavars = cvars;
//        CUDAcommon::handleerror( cudaGetLastError(),"FilamentStretchingHarmonicenergy", "FilamentStretchingHarmonic.cu");
//        double* gpu_Utot = CUDAcommon::getCUDAvars().gpu_energy;
//        addvector<<<1,1,0,stream>>>(gU_i, params, gU_sum, gpu_Utot);
//        CUDAcommon::handleerror( cudaGetLastError() ,"FilamentStretchingHarmonicenergy", "FilamentStretchingHarmonic.cu");
    }

    if(blocksnthreadsez[1]>0) {
        FilamentStretchingHarmonicenergyz << < blocksnthreadsez[0], blocksnthreadsez[1], (12 * blocksnthreadsez[1]) *
                                                sizeof(double), stream>> > (coord, f, beadSet, kstr, eql, params, gU_i,
                                                z, CUDAcommon::getCUDAvars().gculpritID,
                                                CUDAcommon::getCUDAvars().gculpritFF,
                                                CUDAcommon::getCUDAvars().gculpritinteraction, gFF, ginteraction);
        CUDAcommon::handleerror( cudaGetLastError(),"FilamentStretchingHarmonicenergyz", "FilamentStretchingHarmonic"
                ".cu");
//        auto cvars = CUDAcommon::getCUDAvars();
//        cvars.streamvec.push_back(&stream);
//        CUDAcommon::cudavars = cvars;
//        CUDAcommon::handleerror(cudaGetLastError(),"FilamentStretchingHarmonicenergyz", "FilamentStretchingHarmonic"
//                ".cu");
//        double* gpu_Utot = CUDAcommon::getCUDAvars().gpu_energy;
//        addvector<<<1,1,0,stream>>>(gU_i,params, gU_sum, gpu_Utot);
//        CUDAcommon::handleerror(cudaGetLastError(),"FilamentStretchingHarmonicenergyz", "FilamentStretchingHarmonic"
//                ".cu");
    }
    if(blocksnthreadse[1]<=0 && blocksnthreadsez[1]<=0)
        return NULL;
    else{
        auto cvars = CUDAcommon::getCUDAvars();
        cvars.streamvec.push_back(&stream);
        CUDAcommon::cudavars = cvars;
        double* gpu_Utot = CUDAcommon::getCUDAvars().gpu_energy;
//        addvector<<<1,1,0,stream>>>(gU_i, params, gU_sum, gpu_Utot);
//        cudaStreamSynchronize(stream);
//        addvectorred<<<1,200,200*sizeof(double),stream>>>(gU_i,params, gU_sum, gpu_Utot);
//            double Sum[1];
//        CUDAcommon::handleerror(cudaMemcpy(Sum, gU_sum, sizeof(double), cudaMemcpyDeviceToHost));
//        double *gU_sum2;
//        CUDAcommon::handleerror(cudaMalloc((void **) &gU_sum2, sizeof(double)));
//        std::cout<<"bntaddvec "<<bntaddvec2.at(0)<<" "<<bntaddvec2.at(1)<<" "<<bntaddvec2.at(0)<<" "
//                ""<<bntaddvec2.at(2)<<" "<<bntaddvec2.at(3)<<endl;
        resetdoublevariableCUDA<<<1,1,0,stream>>>(gU_sum);
        addvectorred2<<<bntaddvec2.at(2),bntaddvec2.at(3), bntaddvec2.at(3) * sizeof(double),stream>>>(gU_i,
                params, gU_sum, gpu_Utot);
//        double Sum2[1];
//        CUDAcommon::handleerror(cudaMemcpy(Sum2, gU_sum2, sizeof(double), cudaMemcpyDeviceToHost));
//        cudaFree(gU_sum2);
//        std::cout<<Sum[0]<<" "<<Sum2[0]<<endl;
//        cudaStreamSynchronize(stream);
        CUDAcommon::handleerror( cudaGetLastError() ,"FilamentStretchingHarmonicenergy", "FilamentStretchingHarmonic.cu");
        return gU_sum;}
}

void FilamentStretchingHarmonic::forces(double *coord, double *f, int *beadSet,
                                          double *kstr, double *eql, int *params){
    if(blocksnthreadsf[1]>0) {
        FilamentStretchingHarmonicforces << < blocksnthreadsf[0], blocksnthreadsf[1], (6 * blocksnthreadsf[1]) *
                                                sizeof(double), stream >> > (coord, f, beadSet, kstr, eql, params);
        auto cvars = CUDAcommon::getCUDAvars();
        cvars.streamvec.push_back(&stream);
        CUDAcommon::cudavars = cvars;
        CUDAcommon::handleerror(cudaGetLastError(),"FilamentStretchingHarmonicforces", "FilamentStretchingHarmonic.cu");
//        CUDAcommon::handleerror(cudaDeviceSynchronize());
    }
}

 void FilamentStretchingHarmonic::checkforculprit() {
        CUDAcommon::printculprit("FilamentStretching","FilamentStretchingHarmonic");
        Filament* f;
        f = (Filament*)(Cylinder::getCylinders()[CUDAcommon::getCUDAvars().culpritID[0]]->getParent());
        cout<<"Printing culprit Filament information."<<endl;
        f->printSelf();
        exit(EXIT_FAILURE);
}
#endif
double FilamentStretchingHarmonic::energy(double *coord, double *f, int *beadSet,
                                          double *kstr, double *eql){


    int n = FilamentStretching<FilamentStretchingHarmonic>::n;
    int nint = Cylinder::getCylinders().size();

    double *coord1, *coord2, dist, U_i;

    double U = 0.0;

    for(int i = 0; i < nint; i += 1) {

        coord1 = &coord[3 * beadSet[n * i]];
        coord2 = &coord[3 * beadSet[n * i + 1]];
        dist = twoPointDistance(coord1, coord2) - eql[i];

        U_i = 0.5 * kstr[i] * dist * dist;
//        std::cout<<"S "<<i<<" "<<dist<<" "<<U_i<<endl;
        if(fabs(U_i) == numeric_limits<double>::infinity()
           || U_i != U_i || U_i < -1.0) {

            //set culprit and return

            FilamentInteractions::_filamentCulprit = (Filament*)(Cylinder::getCylinders()[i]->getParent());

            return -1;
        }

        U += U_i;
    }
    return U;
}

double FilamentStretchingHarmonic::energy(double *coord, double * f, int *beadSet,
                                          double *kstr, double *eql, double d){

    int n = FilamentStretching<FilamentStretchingHarmonic>::n;
    int nint = Cylinder::getCylinders().size();

    double *coord1, *coord2, *f1, *f2, dist;
    double *v1 = new double[3];
    double *v2 = new double[3];

    double U = 0.0;

    for(int i = 0; i < nint; i += 1) {

        coord1 = &coord[3 * beadSet[n * i]];
        coord2 = &coord[3 * beadSet[n * i + 1]];

        f1 = &f[3 * beadSet[n * i]];
        f2 = &f[3 * beadSet[n * i + 1]];

        dist = twoPointDistanceStretched(coord1, f1,  coord2, f2, d) - eql[i];
//        std::cout<<"S "<<i<<" "<<dist<<" "<<0.5 * kstr[i] * dist * dist<<endl;
        U += 0.5 * kstr[i] * dist * dist;
    }
    delete v1;
    delete v2;

    return U;
}

void FilamentStretchingHarmonic::forces(double *coord, double *f, int *beadSet,
                                        double *kstr, double *eql){


    int n = FilamentStretching<FilamentStretchingHarmonic>::n;
    int nint = Cylinder::getCylinders().size();

    double *coord1, *coord2, dist, invL;
    double f0, *f1, *f2;

    for(int i = 0; i < nint; i += 1) {
        coord1 = &coord[3 * beadSet[n * i]];
        coord2 = &coord[3 * beadSet[n * i + 1]];

        dist = twoPointDistance(coord1, coord2);
        invL = 1 / dist;

        f0 = kstr[i] * ( dist - eql[i] ) * invL;

        f1 = &f[3 * beadSet[n * i]];
        f2 = &f[3 * beadSet[n * i + 1]];

        f2[0] +=  f0 * ( coord1[0] - coord2[0] );
        f2[1] +=  f0 * ( coord1[1] - coord2[1] );
        f2[2] +=  f0 * ( coord1[2] - coord2[2] );

        // force i-1
        f1[0] +=  f0 * ( coord2[0] - coord1[0] );
        f1[1] +=  f0 * ( coord2[1] - coord1[1] );
        f1[2] +=  f0 * ( coord2[2] - coord1[2] );
//        std::cout<<i<<" "<< f0 * ( coord1[0] - coord2[0] )<<" "<<
//        f0 * ( coord1[1] - coord2[1] )<<" "<<
//        f0 * ( coord1[2] - coord2[2] )<<" "<<
//        f0 * ( coord2[0] - coord1[0] )<<" "<<
//        f0 * ( coord2[1] - coord1[1] )<<" "<<
//        f0 * ( coord2[2] - coord1[2] )<<endl;
//        std::cout<<i<<" "<<f1[0]<<" "<<f1[1]<<" "<<f1[2]<<" "<<f2[0]<<" "<<f2[1]<<" "<<f2[2]<<endl;
//        std::cout<<i<<" "<<kstr[i]<<" "<<f0<<" "<<coord1[0]<<" "<<coord1[1]<<" "<<coord1[2]<<" "<<coord2[0]<<" "
//                ""<<coord2[1]<<" "<<coord2[2]<<endl;
    }

}


