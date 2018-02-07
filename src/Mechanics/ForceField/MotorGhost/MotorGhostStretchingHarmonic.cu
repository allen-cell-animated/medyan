
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

#include "MotorGhostStretchingHarmonic.h"
#include "MotorGhostStretching.h"
#include "MotorGhost.h"
#include "MotorGhostStretchingHarmonicCUDA.h"
#include "Bead.h"
#include "MathFunctions.h"
#include "common.h"
#include <cuda.h>
#include <cuda_runtime.h>
#include "nvToolsExt.h"

using namespace mathfunc;
#ifdef CUDAACCL
void MotorGhostStretchingHarmonic::deallocate(){
    CUDAcommon::handleerror(cudaStreamDestroy(stream));
    CUDAcommon::handleerror(cudaFree(gU_i));
    CUDAcommon::handleerror(cudaFree(gU_sum));
    CUDAcommon::handleerror(cudaFree(gFF));
    CUDAcommon::handleerror(cudaFree(ginteraction));
}
void MotorGhostStretchingHarmonic::checkforculprit() {
    CUDAcommon::printculprit("MotorGhostStretching","MotorGhostStretchingHarmonic");
    MotorGhost* m;
    m = MotorGhost::getMotorGhosts()[CUDAcommon::getCUDAvars().culpritID[0]];
    cout<<"Printing culprit Filament information."<<endl;
    m->printSelf();
    exit(EXIT_FAILURE);
}
void MotorGhostStretchingHarmonic::optimalblocksnthreads( int nint){
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
        cudaOccupancyMaxPotentialBlockSizeVariableSMem(&minGridSize, &blockSize,
                                                       MotorGhostStretchingHarmonicenergy, blockToSmem, 0);
//    std::cout<<(nint +blockSize -1) / blockSize<<" "<<blockSize<<endl;
//
//    cudaOccupancyMaxPotentialBlockSize( &minGridSize, &blockSize,
//                                        CUDAExclVolRepulsionenergy, 0, 0);
        blocksnthreadse.push_back((nint + blockSize - 1) / blockSize);
        blocksnthreadse.push_back(blockSize);
//    std::cout<<(nint +blockSize -1) / blockSize<<" "<<blockSize<<endl;
        blockSize = 0;

        cudaOccupancyMaxPotentialBlockSizeVariableSMem(&minGridSize, &blockSize,
                                                       MotorGhostStretchingHarmonicenergyz, blockToSmemez, 0);
        blocksnthreadsez.push_back((nint + blockSize - 1) / blockSize);
        blocksnthreadsez.push_back(blockSize);
        blockSize = 0;

        cudaOccupancyMaxPotentialBlockSizeVariableSMem(&minGridSize, &blockSize,
                                                       MotorGhostStretchingHarmonicforces, blockToSmem, 0);
        blocksnthreadsf.push_back((nint + blockSize - 1) / blockSize);
        blocksnthreadsf.push_back(blockSize);

        CUDAcommon::handleerror(cudaMalloc((void **) &gU_i, nint*sizeof(double)));
	    CUDAcommon::handleerror(cudaMalloc((void **) &gU_sum, sizeof(double)));
        char a[] = "MotorGhostFF";
        char b[] = "MotorGhost Stretching Harmonic";
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
#endif
double* MotorGhostStretchingHarmonic::energy(double *coord, double *f, int *beadSet,
                                            double *kstr, double *eql, double *pos1, double *pos2,
                                            int *params) {
    if(blocksnthreadse[1]>0) {
//        nvtxRangePushA("cmsestream");
//        cudaStream_t  stream;
//        cudaEvent_t event;
//        CUDAcommon::handleerror(cudaStreamCreate(&stream));
//        CUDAcommon::handleerror(cudaEventCreate(&event));
//        auto cvars = CUDAcommon::getCUDAvars();
//        cvars.eventvec.push_back(event);
//        cvars.streamvec.push_back(stream);
//        CUDAcommon::cudavars = cvars;
//        nvtxRangePop();

//    double *gU_ii;
    //double *gU_i;
//    double *gc1, *gc2, *gcheckU;
//    double U_ii[blocksnthreads[0] * blocksnthreads[1]];
//    double c1[3*blocksnthreads[0]*blocksnthreads[1]], c2[3*blocksnthreads[0]*blocksnthreads[1]];
    //double U_i[blocksnthreadse[0]*blocksnthreadse[1]];
//    double checkU[blocksnthreads[1]];

//    double ccoord[3*Bead::getBeads().size()];
//    cudaMemcpy(ccoord, coord, 3*Bead::getBeads().size()*sizeof(double), cudaMemcpyDeviceToHost);
//    double cforce[3*Bead::getBeads().size()];
//    cudaMemcpy(cforce, f, 3*Bead::getBeads().size()*sizeof(double), cudaMemcpyDeviceToHost);
//    for(auto i =0; i < Bead::getBeads().size(); i++)
//        std::cout<<ccoord[3 * i]<<" "<<ccoord[3 * i +1]<<" "<<ccoord[3 * i +2]<<" "<<cforce[3 * i]<<" "<<cforce[3 * i
//                                                                                                                +1]<<" "<<cforce[3 * i +2]<<endl;
//    std::cout<<"C+F---------------------------- "<<endl;

//    std::cout<<"MSE Number of Blocks: "<<blocksnthreads[0]<<endl;
//    std::cout<<"Threads per block: "<<blocksnthreads[1]<<endl;

    //TODO  since the number of threads needed is constant through out the minimization, consider storing the pointer.
    //`CUDAcommon::handleerror(cudaMalloc((void **) &gU_i, blocksnthreadse[0]*blocksnthreadse[1]*sizeof(double)));
//    CUDAcommon::handleerror(cudaMalloc((void **) &gU_ii, blocksnthreads[0]*blocksnthreads[1]*sizeof(double)));
//    CUDAcommon::handleerror(cudaMalloc((void **) &gc1, 3*blocksnthreads[1]*sizeof(double)));
//    CUDAcommon::handleerror(cudaMalloc((void **) &gc2, 3*blocksnthreads[1]*sizeof(double)));
//    CUDAcommon::handleerror(cudaMalloc((void **) &gcheckU, blocksnthreads[1]*sizeof(double)));
    //
//    std::cout<<"MSE CUDA"<<endl;
        nvtxRangePushA("cmse");

    MotorGhostStretchingHarmonicenergy<<<blocksnthreadse[0], blocksnthreadse[1], (12 * blocksnthreadse[1]) * sizeof
                                                                                                                  (double), stream>>>
            (coord, f, beadSet, kstr, eql, pos1, pos2, params, gU_i);
        auto cvars = CUDAcommon::getCUDAvars();
        cvars.streamvec.push_back(&stream);
        CUDAcommon::cudavars = cvars;
//        cudaEventRecord(event, stream);
        nvtxRangePop();
        nvtxRangePushA("cmseError");
                CUDAcommon::handleerror( cudaGetLastError(), "MotorGhostStretchingHarmonicenergy",
                                         "MotorGhostStretchingHarmonic.cu");
        nvtxRangePop();
        //    CUDAcommon::handleerror( cudaPeekAtLastError() );
//    CUDAcommon::handleerror( cudaDeviceSynchronize() );
        nvtxRangePushA("cmseadd");
        double* gpu_Utot = CUDAcommon::getCUDAvars().gpu_energy;
        addvector<<<1,1,0,stream>>>(gU_i,params, gU_sum, gpu_Utot, CUDAcommon::getCUDAvars().gculpritID,
                CUDAcommon::getCUDAvars().gculpritFF,
                CUDAcommon::getCUDAvars().gculpritinteraction, gFF, ginteraction);
        nvtxRangePop();
        nvtxRangePushA("cmseError");
        CUDAcommon::handleerror( cudaGetLastError() , "MotorGhostStretchingHarmonicenergy",
                                 "MotorGhostStretchingHarmonic.cu");
        nvtxRangePop();
//        nvtxRangePushA("cmsecopy");
//        double U_i[1];
//        CUDAcommon::handleerror(cudaMemcpy(U_i, gU_sum, sizeof(double),
//                                           cudaMemcpyDeviceToHost));
//        nvtxRangePop();


//        CUDAcommon::handleerror(cudaMemcpy(U_ii, gU_ii, blocksnthreads[0] * blocksnthreads[1] * sizeof(double),
//                                           cudaMemcpyDeviceToHost));
//        CUDAcommon::handleerror(cudaMemcpy(c1, gc1, 3*blocksnthreads[1]*sizeof(double), cudaMemcpyDeviceToHost));
//        CUDAcommon::handleerror(cudaMemcpy(c2, gc2, 3*blocksnthreads[1]*sizeof(double), cudaMemcpyDeviceToHost));
//        CUDAcommon::handleerror(cudaMemcpy(checkU, gcheckU, blocksnthreads[1]*sizeof(double), cudaMemcpyDeviceToHost));

//    for(auto j=0;j<blocksnthreads[0]* blocksnthreads[1];j++) {
//        std::cout << U_ii[j] <<" "<<checkU[j]<<endl;
//    }

//    for(auto j=0;j<blocksnthreads[0]* blocksnthreads[1];j++) {
//        std::cout << U_i[j] << endl;
//    }
//    CylinderExclVolume<CylinderExclVolRepulsion>::numInteractions
//    for(auto j=0;j<blocksnthreads[0] * blocksnthreads[1];j++) {
////        std::cout << c1[3 * j] << " " << c1[3 * j + 1] << " " << c1[3 * j + 2] << " " << c2[3 * j] << " "
////                  << c2[3 * j + 1] << " " << c2[3 * j + 2] << " ";
//        std::cout << U_i[j] << endl;
//    }
//    std::cout<<"*************"<<endl;

//    auto j=blocksnthreads[1]; int idx=0;
//    while(j!=1)
//    {std::cout<<checkU[idx]<<endl;idx=idx+1;j=j/2;}
//    for (auto i=0;i<blocksnthreads[1];i++)
//        std::cout<<checkU[i]<<endl;
//    std::cout<<"----------------"<<endl;


//    if(U_i[0]!=-1.0) {
//        for (auto i = 1; i < blocksnthreadse[0] * blocksnthreadse[1]; i++) {
////            std::cout<<U_i[i]<<endl;
//            U_i[0] = U_i[0] + U_i[i];
//            if (U_i[i] == -1.0) {
//                U_i[0] = -1.0;
//                break;
//            }
//        }
//    }
//    std::cout<<"MS Total energy CUDA   "<<U_i[0]<<endl;

    //CUDAcommon::handleerror(cudaFree(gU_i));

//        CUDAcommon::handleerror(cudaFree(gU_ii));
//        CUDAcommon::handleerror(cudaFree(gc1));
//        CUDAcommon::handleerror(cudaFree(gc2));
//        CUDAcommon::handleerror(cudaFree(gcheckU));
//    gU_i = NULL;
//    gU_ii = NULL;
//    gc1 = NULL;
//    gc2 = NULL;
//    gcheckU = NULL;
//    free(U_i);


//        CUDAcommon::handleerror(cudaStreamSynchronize(stream));
    return gU_sum;}
    else
        return NULL;
}


double* MotorGhostStretchingHarmonic::energy(double *coord, double *f, int *beadSet,
                                            double *kstr, double *eql, double *pos1, double *pos2, double *z,
                                            int *params) {

    ///TEST CODE ///
//    double *gcheckvar;
//    double checkvar[blocksnthreads[0]*blocksnthreads[1]];
//    CUDAcommon::handleerror(cudaMalloc((void **) &gcheckvar, blocksnthreads[0] * blocksnthreads[1]*sizeof(double)));
//    testifitworks<<<blocksnthreads[0], blocksnthreads[1] >>>(coord,gcheckvar);
//    cudaMemcpy(checkvar, gcheckvar, blocksnthreads[0] * blocksnthreads[1]*sizeof(double), cudaMemcpyDeviceToHost);
//    for (auto i = 1; i < blocksnthreads[0] * blocksnthreads[1]; i++) {
//        std::cout<<checkvar[i]<<" ";
//    }
//    std::cout<<endl;
//    cudaFree(gcheckvar);
//    gcheckvar = NULL;
    ///@@END@@///
    if(blocksnthreadsez[1]>0) {
//        cudaStream_t  stream;
//        cudaEvent_t event;
//        CUDAcommon::handleerror(cudaStreamCreate(&stream));
//        CUDAcommon::handleerror(cudaEventCreate(&event));
//        auto cvars = CUDAcommon::getCUDAvars();
//        cvars.eventvec.push_back(event);
//        cvars.streamvec.push_back(stream);
//        CUDAcommon::cudavars = cvars;
//        double dd[1];
//        CUDAcommon::handleerror(cudaMemcpy(dd, z, sizeof(double), cudaMemcpyDeviceToHost));
//        std::cout << "d = " << dd[0] << endl;

//        double *gU_ii;
        //double *gU_i;
//        double *gc1, *gc2, *gcheckU;
//        double U_ii[blocksnthreads[0] * blocksnthreads[1]];
//        double c1[3 * blocksnthreads[0] * blocksnthreads[1]], c2[3 * blocksnthreads[0] * blocksnthreads[1]];
        //double U_i[blocksnthreadsez[0] * blocksnthreadsez[1]];
//        double checkU[33 * blocksnthreads[0] * blocksnthreads[1]];



//    double ccoord[3*Bead::getBeads().size()];
//    cudaMemcpy(ccoord, coord, 3*Bead::getBeads().size()*sizeof(double), cudaMemcpyDeviceToHost);
//    double cforce[3*Bead::getBeads().size()];
//    cudaMemcpy(cforce, f, 3*Bead::getBeads().size()*sizeof(double), cudaMemcpyDeviceToHost);
//        int cparams[2];
//        CUDAcommon::handleerror(cudaMemcpy(cparams, params, 2 * sizeof(int), cudaMemcpyDeviceToHost));
//    for(auto i =0; i < Bead::getBeads().size(); i++)
//
//        std::cout<<ccoord[3 * i]<<" "<<ccoord[3 * i +1]<<" "<<ccoord[3 * i +2]<<" "<<cforce[3 * i]<<" "<<cforce[3 * i
//                                                                                                                +1]<<" "<<cforce[3 * i +2]<<endl;
//        std::cout << cparams[0] << " " << cparams[1] << endl;

//    std::cout<<"C+F Z---------------------------- "<<endl;

//        std::cout << "MSEZ Number of Blocks: " << blocksnthreads[0] << endl;
//        std::cout << "Threads per block: " << blocksnthreads[1] << endl;



        //TODO  since the number of threads needed is constant through out the minimization, consider storing the pointer.
        //CUDAcommon::handleerror(cudaMalloc((void **) &gU_i, blocksnthreadsez[0] * blocksnthreadsez[1] * sizeof (double)));
//        CUDAcommon::handleerror(cudaMalloc((void **) &gU_ii, blocksnthreads[0] * blocksnthreads[1] * sizeof(double)));
//        CUDAcommon::handleerror(cudaMalloc((void **) &gc1, 3 * blocksnthreads[1] * sizeof(double)));
//        CUDAcommon::handleerror(cudaMalloc((void **) &gc2, 3 * blocksnthreads[1] * sizeof(double)));
//        CUDAcommon::handleerror(cudaMalloc((void **) &gcheckU, 33 * blocksnthreads[0] * blocksnthreads[1] * sizeof
//        (double)));
        //
//    size_t freeMem, totalMem;
//
//    cudaMemGetInfo(&freeMem, &totalMem);
//
//    std::cout<<"Memory "<<freeMem<<" "<<totalMem<<endl;
//        struct cudaDeviceProp properties;
//        cudaGetDeviceProperties(&properties, 0);
//        cout << "using " << properties.multiProcessorCount << " multiprocessors" << endl;
//        cout << "max threads per processor: " << properties.maxThreadsPerMultiProcessor << endl;
//        std::cout << 24 * blocksnthreads[1] * sizeof(double) << endl;

//        if(blocksnthreadsez[1]==THREADSPERBLOCK) {
//            blocksnthreadsez[0] = 2 * blocksnthreadsez[0];
//            blocksnthreadsez[1] = 1/2 * blocksnthreadsez[1];
//        }
        nvtxRangePushA("cmsez");
        MotorGhostStretchingHarmonicenergyz << < blocksnthreadsez[0], blocksnthreadsez[1], (24 * blocksnthreadsez[1]) *
                                                                                                sizeof
                (double), stream>> > (coord, f, beadSet, kstr, eql, pos1, pos2, params, gU_i, z );
        auto cvars = CUDAcommon::getCUDAvars();
        cvars.streamvec.push_back(&stream);
        CUDAcommon::cudavars = cvars;
        nvtxRangePop();
        nvtxRangePushA("cmsezError");
        CUDAcommon::handleerror(cudaGetLastError(), "MotorGhostStretchingHarmonicenergyz",
                                "MotorGhostStretchingHarmonic.cu");
        nvtxRangePop();
//    CUDAcommon::handleerror( cudaPeekAtLastError() );
//        CUDAcommon::handleerror(cudaDeviceSynchronize());
        nvtxRangePushA("cmsezadd");
        double* gpu_Utot = CUDAcommon::getCUDAvars().gpu_energy;

        addvector<<<1,1,0,stream>>>(gU_i,params, gU_sum, gpu_Utot, CUDAcommon::getCUDAvars().gculpritID,
                CUDAcommon::getCUDAvars().gculpritFF,
                CUDAcommon::getCUDAvars().gculpritinteraction, gFF, ginteraction);
//        cudaEventRecord(event, stream);
        nvtxRangePop();
        nvtxRangePushA("cmsezError");
        CUDAcommon::handleerror(cudaGetLastError(), "MotorGhostStretchingHarmonicenergyz",
                                "MotorGhostStretchingHarmonic.cu");
        nvtxRangePop();

//        nvtxRangePushA("cmsezcopy");
//        double U_i[1];
//        CUDAcommon::handleerror(cudaMemcpy(U_i, gU_sum, sizeof(double),
//                                           cudaMemcpyDeviceToHost));
//        nvtxRangePop();
//
//        CUDAcommon::handleerror(cudaMemcpy(U_i, gU_i, blocksnthreadsez[0] * blocksnthreadsez[1] * sizeof(double),
//                                                                                                       cudaMemcpyDeviceToHost));
//        CUDAcommon::handleerror(cudaMemcpy(U_ii, gU_ii, blocksnthreads[0] * blocksnthreads[1] * sizeof(double),
//                                           cudaMemcpyDeviceToHost));
//        CUDAcommon::handleerror(cudaMemcpy(c1, gc1, 3 * blocksnthreads[1] * sizeof(double), cudaMemcpyDeviceToHost));
//        CUDAcommon::handleerror(cudaMemcpy(c2, gc2, 3 * blocksnthreads[1] * sizeof(double), cudaMemcpyDeviceToHost));
//        CUDAcommon::handleerror(cudaMemcpy(checkU, gcheckU, 33 * blocksnthreads[0] * blocksnthreads[1] * sizeof(double),
//                                           cudaMemcpyDeviceToHost));
//
//        for(auto i=0;i<blocksnthreads[0] * blocksnthreads[1]; i++) {
//            for (auto iter = 0; iter < 33; iter++)
//                std::cout << checkU[33 * i + iter] << " ";
//            std::cout<<U_i[i]<<endl;
//        }
//        if (U_i[0] != -1.0) {
//            for (auto i = 1; i < blocksnthreadsez[0] * blocksnthreadsez[1]; i++) {
//                U_i[0] = U_i[0] + U_i[i];
//                if (U_i[i] == -1.0) {
//                    U_i[0] = -1.0;
//                    break;
//                }
//            }
//        }
//        std::cout << "MSZ Total energy CUDA   " << U_i[0] << endl;

//        CUDAcommon::handleerror(cudaFree(gU_i));

//        CUDAcommon::handleerror(cudaFree(gU_ii));
//        CUDAcommon::handleerror(cudaFree(gc1));
//        CUDAcommon::handleerror(cudaFree(gc2));
//        CUDAcommon::handleerror(cudaFree(gcheckU));
//    gU_i = NULL;
//    gU_ii = NULL;
//    gc1 = NULL;
//    gc2 = NULL;
//    gcheckU = NULL;
//    free(U_i);


//        CUDAcommon::handleerror(cudaStreamSynchronize(stream));
        return gU_sum;
    }else
        return NULL;
}

double MotorGhostStretchingHarmonic::energy(double *coord, double *f, int *beadSet,
                                            double *kstr, double *eql, double *pos1, double *pos2) {

    int n = MotorGhostStretching<MotorGhostStretchingHarmonic>::n;
    int nint = MotorGhost::getMotorGhosts().size();

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
            MotorGhostInteractions::_motorCulprit = MotorGhost::getMotorGhosts()[i];

            return -1;
        }

        U += U_i;
//        std::cout<<U_i<<endl;
    }
//    std::cout<<"MS Total energy serial "<< U <<endl;
    delete v1;
    delete v2;

    return U;
}

double MotorGhostStretchingHarmonic::energy(double *coord, double * f, int *beadSet,
                                            double *kstr, double *eql, double *pos1, double *pos2, double d){

    int n = MotorGhostStretching<MotorGhostStretchingHarmonic>::n;
    int nint = MotorGhost::getMotorGhosts().size();

    double *coord1, *coord2, *coord3, *coord4, *f1, *f2, *f3, *f4, dist, U_i;
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
        U_i = 0.5 * kstr[i] * dist * dist;

//        std::cout<<v1[0]<<" "<<v1[1]<<" "<<v1[2]<<" "<<v2[0]<<" "<<v2[1]<<" "<<v2[2]<<" "<<coord1[0]<<" "
//                ""<<coord1[1]<<" "<<coord1[2]<<" "<<coord3[0]<<" "<<coord3[1]<<" "<<coord3[2]<<" "<<f1[0]<<" "
//                ""<<f1[1]<<" "<<f1[2]<<" "<<f3[0]<<" "<<f3[1]<<" "<<f3[2]<<" "<<coord2[0]<<" "<<coord2[1]<<" "
//                ""<<coord2[2]<<" "<<coord4[0]<<" "<<coord4[1]<<" "<<coord4[2]<<" "<<f2[0]<<" "<<f2[1]<<" "<<f2[2]<<" "
//                ""<<f4[0]<<" "<<f4[1]<<" "<<f4[2]<<" "<<pos1[i]<<" "<<pos2[i]<<" "<<d<<" "<<U_i<<endl;
        if(fabs(U_i) == numeric_limits<double>::infinity()
           || U_i != U_i || U_i < -1.0) {

            //set culprit and return
            MotorGhostInteractions::_motorCulprit = MotorGhost::getMotorGhosts()[i];

            return -1;
        }

        U += U_i;
    }

//    std::cout<<"MS Total energy serial "<< U <<endl;
    delete v1;
    delete v2;

    return U;

}

void MotorGhostStretchingHarmonic::forces(double *coord, double *f, int *beadSet,
                                          double *kstr, double *eql, double *pos1, double *pos2, int *params){
//    cudaEvent_t start, stop;
//    CUDAcommon::handleerror(cudaEventCreate( &start));
//    CUDAcommon::handleerror(cudaEventCreate( &stop));
//    CUDAcommon::handleerror(cudaEventRecord( start, 0));
    if(blocksnthreadsf[1]>0) {
//        double *gU_i;
//        double *gc1, *gc2, *gcheckU;

//        double c1[3 * blocksnthreads[0] * blocksnthreads[1]], c2[3 * blocksnthreads[0] * blocksnthreads[1]];
//        double cvar[36 * blocksnthreads[0] * blocksnthreads[1]];

//        std::cout << "MSF Number of Blocks: " << blocksnthreads[0] << endl;
//        std::cout << "Threads per block: " << blocksnthreads[1] << endl;

        //TODO  since the number of threads needed is constant through out the minimization, consider storing the pointer.
//        CUDAcommon::handleerror(cudaMalloc((void **) &gU_i, 36 * blocksnthreads[0] * blocksnthreads[1] * sizeof
//                                                                                                                (double)));


//        double F_c[3*Bead::getBeads().size()];
//        double C_c[3*Bead::getBeads().size()];
//        //TODO remove this later need not copy forces back to CPU.
//        CUDAcommon::handleerror(cudaMemcpy(F_c, f, 3 * Bead::getBeads().size() *sizeof(double),
//                                           cudaMemcpyDeviceToHost));
//        CUDAcommon::handleerror(cudaMemcpy(C_c, coord, 3 * Bead::getBeads().size() *sizeof(double),
//                                           cudaMemcpyDeviceToHost));
//        for(int iter=0;iter<Bead::getBeads().size();iter++) {
//            std::cout << C_c[3 * iter] << " " << C_c[3 * iter + 1] << " " << C_c[3 * iter + 2]<<" "<<F_c[3 * iter] <<
//            " " << F_c[3 * iter + 1] << " " << F_c[3 * iter + 2] <<endl;
//        }
//
//        std::cout<<"check ends "<<blocksnthreads[0]<<" "<<blocksnthreads[1]<<endl;
//
//    size_t freeMem, totalMem;
//
//    cudaMemGetInfo(&freeMem, &totalMem);
//
//    std::cout<<"Memory "<<freeMem<<" "<<totalMem<<endl;
//        struct cudaDeviceProp properties;
//        cudaGetDeviceProperties(&properties, 0);
//        cout << "using " << properties.multiProcessorCount << " multiprocessors" << endl;
//        cout << "max threads per processor: " << properties.maxThreadsPerMultiProcessor << endl;
//        std::cout << 36 *  blocksnthreads[0] *blocksnthreads[1] * sizeof(double) << endl;

        MotorGhostStretchingHarmonicforces << < blocksnthreadsf[0], blocksnthreadsf[1], (12 *
        blocksnthreadsf[1]) * sizeof (double), stream >> > (coord, f, beadSet, kstr, eql, pos1, pos2, params);
        auto cvars = CUDAcommon::getCUDAvars();
        cvars.streamvec.push_back(&stream);
        CUDAcommon::cudavars = cvars;
        //CUDAcommon::handleerror(cudaDeviceSynchronize());
        CUDAcommon::handleerror(cudaGetLastError(), "MotorGhostStretchingHarmonicforces",
                                "MotorGhostStretchingHarmonic.cu");

//    CUDAcommon::handleerror( cudaPeekAtLastError() );
//        CUDAcommon::handleerror(cudaDeviceSynchronize());

//        CUDAcommon::handleerror(cudaMemcpy(cvar, gU_i, 36 * blocksnthreads[0]*blocksnthreads[1]*sizeof(double),
//                                           cudaMemcpyDeviceToHost));
//        for(auto i=0; i<blocksnthreads[0]*blocksnthreads[1]; i++) {
//            for(auto iter=0;iter<36;iter++) {
//                std::cout <<cvar[36 * i + iter]<<" ";
//            }
//            std::cout<<endl;
//        }

//    cudaMemcpy(F_i, f, 3 * blocksnthreads[0]*blocksnthreads[1]*sizeof(double), cudaMemcpyDeviceToHost);
//        CUDAcommon::handleerror(cudaFree(gU_i));
    }
//    CUDAcommon::handleerror(cudaEventRecord( stop, 0));
//    CUDAcommon::handleerror(cudaEventSynchronize(stop));
//    float elapsedtime;
//    CUDAcommon::handleerror(cudaEventElapsedTime(&elapsedtime, start, stop));
//    CUDAvars cvars=CUDAcommon::getCUDAvars();
//    cvars.Ccforce += elapsedtime;
//    std::cout<<"C CFM "<<elapsedtime<<endl;
//    CUDAcommon::cudavars=cvars;
//    CUDAcommon::handleerror(cudaEventDestroy(start));
//    CUDAcommon::handleerror(cudaEventDestroy(stop));
}

void MotorGhostStretchingHarmonic::forces(double *coord, double *f, int *beadSet,
                                          double *kstr, double *eql, double *pos1, double *pos2){

    int n = MotorGhostStretching<MotorGhostStretchingHarmonic>::n;
    int nint = MotorGhost::getMotorGhosts().size();

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

//        std::cout<<coord1[0]<<" "<<coord1[1]<<" "<<coord1[2]<<" "<<coord2[0]<<" "<<coord2[1]<<" "<<coord2[2]<<" "
//                ""<<coord3[0]<<" "<<coord3[1]<<" "<<coord3[2]<<" "<<coord4[0]<<" "<<coord4[1]<<" "<<coord4[2]<<" ";
//        std::cout<<-f0 * ( v1[0] - v2[0] ) * (1 - pos1[i])<<" "
//                 <<-f0 * ( v1[1] - v2[1] ) * (1 - pos1[i])<<" "
//                 <<-f0 * ( v1[2] - v2[2] ) * (1 - pos1[i])<<" "
//                 <<-f0 * ( v1[0] - v2[0] ) * (pos1[i])<<" "
//                 <<-f0 * ( v1[1] - v2[1] ) * (pos1[i])<<" "
//                 <<-f0 * ( v1[2] - v2[2] ) * (pos1[i])<<" "
//                 <<f0 * ( v1[0] - v2[0] ) * (1 - pos2[i])<<" "
//                 <<f0 * ( v1[1] - v2[1] ) * (1 - pos2[i])<<" "
//                 <<f0 * ( v1[2] - v2[2] ) * (1 - pos2[i])<<" "
//                 <<f0 * ( v1[0] - v2[0] ) * (pos2[i])<<" "
//                 <<f0 * ( v1[1] - v2[1] ) * (pos2[i])<<" "
//                 <<f0 * ( v1[2] - v2[2] ) * (pos2[i])<<" ";
//        std::cout<<f1[0]<<" "<<f1[1]<<" "<<f1[2]<<" "<<f2[0]<<" "<<f2[1]<<" "<<f2[2]<<" "<<f3[0]<<" "
//                ""<<f3[1]<<" "<<f3[2]<<" "<<f4[0]<<" "<<f4[1]<<" "<<f4[2]<<endl;

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
    }
    delete v1;
    delete v2;

}
