
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v3.2.1
//
//  Copyright (2015-2018)  Papoian Lab, University of Maryland
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

#ifdef CUDAACCL
#include "nvToolsExt.h"
#include <cuda.h>
#include <cuda_runtime.h>
#endif

using namespace mathfunc;
#ifdef CUDAACCL
void FilamentStretchingHarmonic::deallocate(){
//    if(!(CUDAcommon::getCUDAvars().conservestreams))
//        CUDAcommon::handleerror(cudaStreamDestroy(stream));
    CUDAcommon::handleerror(cudaFree(gU_i));
    CUDAcommon::handleerror(cudaFree(gU_sum));
    CUDAcommon::handleerror(cudaFree(gFF));
    CUDAcommon::handleerror(cudaFree(ginteraction));
}
void FilamentStretchingHarmonic::optimalblocksnthreads(int nint, cudaStream_t stream_pass){
    //CUDA stream create
//    if(stream == NULL || !(CUDAcommon::getCUDAvars().conservestreams))
//        CUDAcommon::handleerror(cudaStreamCreate(&stream));
    stream = stream_pass;
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
                                                       FilamentStretchingHarmonicenergyz, blockToSmemZero, 0);
//                                                       FilamentStretchingHarmonicenergyz, blockToSmem, 0);
        blocksnthreadsez.push_back((nint + blockSize - 1) / blockSize);
        blocksnthreadsez.push_back(blockSize);
        blockSize = 0;

        cudaOccupancyMaxPotentialBlockSizeVariableSMem(&minGridSize, &blockSize,
                                                       FilamentStretchingHarmonicforces, blockToSmemF, 0);
//                                                       FilamentStretchingHarmonicforces, blockToSmemZero, 0);
        blocksnthreadsf.push_back((nint + blockSize - 1) / blockSize);
        blocksnthreadsf.push_back(blockSize);
        //get addition vars
        bntaddvec2.clear();
        bntaddvec2 = getaddred2bnt(nint);
        CUDAcommon::handleerror(cudaMalloc((void **) &gU_i, bntaddvec2.at(0)*sizeof(floatingpoint)));
        CUDAcommon::handleerror(cudaMemsetAsync(gU_i, 0, bntaddvec2.at(0) * sizeof(floatingpoint),stream));
        CUDAcommon::handleerror(cudaMalloc((void **) &gU_sum, sizeof(floatingpoint)));
        char a[] = "FilamentFF";
        char b[] = "Filament Stretching Harmonic";
        CUDAcommon::handleerror(cudaMalloc((void **) &gFF, 100 * sizeof(char)));
        CUDAcommon::handleerror(cudaMalloc((void **) &ginteraction, 100 * sizeof(char)));
        CUDAcommon::handleerror(cudaMemcpyAsync(gFF, a, 100 * sizeof(char),
                                            cudaMemcpyHostToDevice, stream));
        CUDAcommon::handleerror(cudaMemcpyAsync(ginteraction, b, 100 * sizeof(char),
                                           cudaMemcpyHostToDevice, stream));
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
floatingpoint* FilamentStretchingHarmonic::energy(floatingpoint *coord, floatingpoint *f, int *beadSet,
                                             floatingpoint *kstr, floatingpoint *eql, int *params) {
//    if(blocksnthreadse[1]>0) {
//        FilamentStretchingHarmonicenergy<<<blocksnthreadse[0], blocksnthreadse[1], (6 * blocksnthreadse[1]) * sizeof
//                (floatingpoint), stream>>> (coord, f, beadSet, kstr, eql, params, gU_i, CUDAcommon::getCUDAvars()
//                .gculpritID,
//                CUDAcommon::getCUDAvars().gculpritFF,
//                CUDAcommon::getCUDAvars().gculpritinteraction, gFF, ginteraction);
//        auto cvars = CUDAcommon::getCUDAvars();
//        cvars.streamvec.push_back(&stream);
//        CUDAcommon::cudavars = cvars;
//        CUDAcommon::handleerror( cudaGetLastError(),"FilamentStretchingHarmonicenergy", "FilamentStretchingHarmonic.cu");
//        floatingpoint* gpu_Utot = CUDAcommon::getCUDAvars().gpu_energy;
//        addvector<<<1,1,0,stream>>>(gU_i, params, gU_sum, gpu_Utot);
//        CUDAcommon::handleerror( cudaGetLastError() ,"FilamentStretchingHarmonicenergy", "FilamentStretchingHarmonic.cu");
//
//        return gU_sum;
//    }
//    else
//        return NULL;
}


floatingpoint* FilamentStretchingHarmonic::energy(floatingpoint *coord, floatingpoint *f, int *beadSet,
                                             floatingpoint *kstr, floatingpoint *eql, floatingpoint *z, int *params) {

//    if(blocksnthreadse[1]>0) {
//        FilamentStretchingHarmonicenergy<<<blocksnthreadse[0], blocksnthreadse[1], (6 * blocksnthreadse[1]) * sizeof
//                (floatingpoint), stream>>> (coord, f, beadSet, kstr, eql, params, gU_i, z, CUDAcommon::getCUDAvars()
//                .gculpritID,
//                CUDAcommon::getCUDAvars().gculpritFF,
//                CUDAcommon::getCUDAvars().gculpritinteraction, gFF, ginteraction);
//        CUDAcommon::handleerror( cudaGetLastError(),"FilamentStretchingHarmonicenergy", "FilamentStretchingHarmonic.cu");
//    }

    if(blocksnthreadsez[1]>0) {
        auto boolvarvec = CUDAcommon::cudavars.backtrackbools;
//        FilamentStretchingHarmonicenergyz << < blocksnthreadsez[0], blocksnthreadsez[1], (12 * blocksnthreadsez[1]) *
        FilamentStretchingHarmonicenergyz << < blocksnthreadsez[0], blocksnthreadsez[1], (0 )*
                                                sizeof(floatingpoint), stream>> > (coord, f, beadSet,
                                                kstr, eql, params, gU_i,
                                                CUDAcommon::cudavars.gpu_energyvec,
                                                z, CUDAcommon::getCUDAvars().gculpritID,
                                                CUDAcommon::getCUDAvars().gculpritFF,
                                                CUDAcommon::getCUDAvars().gculpritinteraction, gFF, ginteraction, boolvarvec.at(0),
                                                boolvarvec.at(1));
    }
    if(blocksnthreadse[1]<=0 && blocksnthreadsez[1]<=0)
        return NULL;
    else{
        auto cvars = CUDAcommon::getCUDAvars();
        cvars.streamvec.push_back(&stream);
        CUDAcommon::cudavars = cvars;
#ifdef CUDA_INDIVIDUAL_ESUM
        floatingpoint* gpu_Utot = CUDAcommon::getCUDAvars().gpu_energy;
        resetfloatingpointvariableCUDA<<<1,1,0,stream>>>(gU_sum);
        addvectorred2<<<bntaddvec2.at(2),bntaddvec2.at(3), bntaddvec2.at(3) * sizeof(floatingpoint),stream>>>(gU_i,
                params, gU_sum, gpu_Utot);
#endif
        CUDAcommon::handleerror( cudaGetLastError() ,"FilamentStretchingHarmonicenergy", "FilamentStretchingHarmonic.cu");
        return gU_sum;}
}

void FilamentStretchingHarmonic::forces(floatingpoint *coord, floatingpoint *f, int *beadSet,
                                          floatingpoint *kstr, floatingpoint *eql, int *params){
    if(blocksnthreadsf[1]>0) {
        FilamentStretchingHarmonicforces << < blocksnthreadsf[0], blocksnthreadsf[1], (6 * blocksnthreadsf[1]) *
//        FilamentStretchingHarmonicforces << < blocksnthreadsf[0], blocksnthreadsf[1], (0) *
                                                sizeof(floatingpoint), stream >> > (coord, f, beadSet, kstr, eql, params);
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
floatingpoint FilamentStretchingHarmonic::energy(floatingpoint *coord, floatingpoint *f, int *beadSet,
                                          floatingpoint *kstr, floatingpoint *eql){


    int n = FilamentStretching<FilamentStretchingHarmonic>::n;
    int nint = Cylinder::getCylinders().size();

    floatingpoint *coord1, *coord2, dist;

    floatingpoint U_i, U = 0.0;

    for(int i = 0; i < nint; i += 1) {

        coord1 = &coord[3 * beadSet[n * i]];
        coord2 = &coord[3 * beadSet[n * i + 1]];
        dist = twoPointDistance(coord1, coord2) - eql[i];

        U_i = 0.5 * kstr[i] * dist * dist;
//        std::cout<<"S "<<i<<" "<<dist<<" "<<U_i<<endl;
        if(fabs(U_i) == numeric_limits<floatingpoint>::infinity()
           || U_i != U_i || U_i < -1.0) {

            //set culprit and return

            FilamentInteractions::_filamentCulprit = (Filament*)(Cylinder::getCylinders()[i]->getParent());

            return -1;
        }

        U += U_i;
    }
    return U;
}

floatingpoint FilamentStretchingHarmonic::energy(floatingpoint *coord, floatingpoint * f, int *beadSet,
                                          floatingpoint *kstr, floatingpoint *eql, floatingpoint d){

    int n = FilamentStretching<FilamentStretchingHarmonic>::n;
    int nint = Cylinder::getCylinders().size();

    floatingpoint *coord1, *coord2, dist;
    floatingpoint *f1, *f2;
    floatingpoint *v1 = new floatingpoint[3];
    floatingpoint *v2 = new floatingpoint[3];

    floatingpoint U = 0.0;

    for(int i = 0; i < nint; i += 1) {

        coord1 = &coord[3 * beadSet[n * i]];
        coord2 = &coord[3 * beadSet[n * i + 1]];

        f1 = &f[3 * beadSet[n * i]];
        f2 = &f[3 * beadSet[n * i + 1]];

        dist = twoPointDistanceStretched(coord1, f1,  coord2, f2, d) - eql[i];
//        std::cout<<"S "<<i<<" "<<dist<<" "<<0.5 * kstr[i] * dist * dist<<endl;
        U += 0.5 * kstr[i] * dist * dist;
    }
    delete [] v1;
    delete [] v2;

    return U;
}

void FilamentStretchingHarmonic::forces(floatingpoint *coord, floatingpoint *f, int *beadSet,
                                        floatingpoint *kstr, floatingpoint *eql){


    int n = FilamentStretching<FilamentStretchingHarmonic>::n;
    int nint = Cylinder::getCylinders().size();

    floatingpoint *coord1, *coord2, dist, invL;
    floatingpoint f0, *f1, *f2;

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

        #ifdef CHECKFORCES_INF_NAN
        if(checkNaN_INF(f1, 0, 2)||checkNaN_INF(f2,0,2)){
            cout<<"Filament Stretching Force becomes infinite. Printing data "<<endl;

            auto cyl = Cylinder::getCylinders()[i];
            cout<<"Cylinder ID "<<cyl->getId()<<" with cindex "<<cyl->getStableIndex()<<
            " and bIndex "<< cyl->getFirstBead()->getIndex()<<" "<<cyl->getSecondBead()
            ->getIndex()<<endl;

            cout<<"Printing coords"<<endl;
            cout<<coord1[0]<<" "<<coord1[1]<<" "<<coord1[2]<<endl;
            cout<<coord2[0]<<" "<<coord2[1]<<" "<<coord2[2]<<endl;
            cout<<"Printing force"<<endl;
            cout<<f1[0]<<" "<<f1[1]<<" "<<f1[2]<<endl;
            cout<<f2[0]<<" "<<f2[1]<<" "<<f2[2]<<endl;
	        cout<<"Printing binary Coords"<<endl;
	        printvariablebinary(coord1,0,2);
	        printvariablebinary(coord2,0,2);
	        cout<<"Printing binary Force"<<endl;
	        printvariablebinary(f1,0,2);
	        printvariablebinary(f2,0,2);
            exit(EXIT_FAILURE);
        }
		#endif
    }
}