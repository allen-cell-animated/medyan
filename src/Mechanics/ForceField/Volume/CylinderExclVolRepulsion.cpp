
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

/*
 Cylinder 1: x1------x2
 Cylinder 2: y1------y2


 x = x2-x1
 y = y2-y1
 z = x1-y1

 a = x.x;
 b = y.y;
 c = z.z;
 d = x.y;
 e = x.z;
 f = y.z;

 */

#include "CylinderExclVolRepulsion.h"
#include "CylinderExclVolRepulsionCUDA.h"
#include "CylinderExclVolume.h"

#include "Bead.h"
#include "Cylinder.h"

#include "MathFunctions.h"
#include "SysParams.h"
#include <limits>
#include "CGMethod.h"
#ifdef CUDAACCL
#include "nvToolsExt.h"
#include <cuda.h>
#include <cuda_runtime.h>
#endif
typedef std::numeric_limits< double > dbl;

using namespace mathfunc;
#ifdef CUDAACCL
//struct unaryfn: std::unary_function<size_t, unsigned long> {
//    int operator()(unsigned long i) const { return 12* i *sizeof(double); }
//
//};
void CylinderExclVolRepulsion::deallocate(){
//    if(!(CUDAcommon::getCUDAvars().conservestreams))
//        CUDAcommon::handleerror(cudaStreamDestroy(stream),"cuda stream",
// "CylinderExclVolumeRepulsion.cu");
    CUDAcommon::handleerror(cudaFree(gU_i),"cudaFree", "CylinderExclVolume.cu");
    CUDAcommon::handleerror(cudaFree(gU_sum),"cudaFree", "CylinderExclVolume.cu");
    CUDAcommon::handleerror(cudaFree(gFF),"cudaFree", "CylinderExclVolume.cu");
    CUDAcommon::handleerror(cudaFree(ginteraction),"cudaFree", "CylinderExclVolume.cu");
}
void CylinderExclVolRepulsion::optimalblocksnthreads( int nint, cudaStream_t stream_pass) {
//    //CUDA stream create
//    if(stream == NULL || !(CUDAcommon::getCUDAvars().conservestreams))
//        CUDAcommon::handleerror(cudaStreamCreate(&stream),"cuda stream", "CylinderExclVolume.cu");
    //
    stream = stream_pass;
    blocksnthreadse.clear();
    blocksnthreadsez.clear();
    blocksnthreadsf.clear();
    if(nint>0){
    int blockSize;   // The launch configurator returned block size
    int minGridSize; // The minimum grid size needed to achieve the
    // maximum occupancy for a full device launch
    int gridSize;    // The actual grid size needed, based on input size
//    unaryfn::argument_type blksize;
//    unaryfn::result_type result;
//    unaryfn ufn;

    CUDAcommon::handleerror(cudaOccupancyMaxPotentialBlockSizeVariableSMem(&minGridSize, &blockSize,
                            CUDAExclVolRepulsionenergy, blockToSmem, 0),"cuda occupancy", "CylinderExclVolume.cu");
//    std::cout<<(nint +blockSize -1) / blockSize<<" "<<blockSize<<endl;
//
//    cudaOccupancyMaxPotentialBlockSize( &minGridSize, &blockSize,
//                                        CUDAExclVolRepulsionenergy, 0, 0);
    blocksnthreadse.push_back((nint + blockSize - 1) / blockSize);
    blocksnthreadse.push_back(blockSize);
//    std::cout<<(nint +blockSize -1) / blockSize<<" "<<blockSize<<endl;
    blockSize = 0;

        CUDAcommon::handleerror(cudaOccupancyMaxPotentialBlockSizeVariableSMem(&minGridSize, &blockSize,
                                CUDAExclVolRepulsionenergyz, blockToSmem, 0),"cuda occupancy", "CylinderExclVolume.cu");
    blocksnthreadsez.push_back((nint + blockSize - 1) / blockSize);
    blocksnthreadsez.push_back(blockSize);
    blockSize = 0;

        CUDAcommon::handleerror(cudaOccupancyMaxPotentialBlockSizeVariableSMem(&minGridSize, &blockSize,
                                CUDAExclVolRepulsionforce, blockToSmem, 0),"cuda occupancy", "CylinderExclVolume.cu");
    blocksnthreadsf.push_back((nint + blockSize - 1) / blockSize);
    blocksnthreadsf.push_back(blockSize);
//get addition vars
        bntaddvec2.clear();
        bntaddvec2 = getaddred2bnt(nint);
        CUDAcommon::handleerror(cudaMalloc((void **) &gU_i, bntaddvec2.at(0)*sizeof(double)));
        CUDAcommon::handleerror(cudaMemsetAsync(gU_i, 0, bntaddvec2.at(0) * sizeof(double), stream),
                                "cuda data transfer",
                                "CylinderExclVolume.cu");
        CUDAcommon::handleerror(cudaMalloc((void **) &gU_i, nint*sizeof(double)),"cuda data transfer",
                                "CylinderExclVolume.cu");
        CUDAcommon::handleerror(cudaMalloc((void **) &gU_sum, sizeof(double)),"cuda data transfer",
                                "CylinderExclVolume.cu");
        char a[] = "Excluded Volume";
        char b[] =  "Cylinder Excluded Volume";
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
double* CylinderExclVolRepulsion::energy(double *coord, double *f, int *beadSet,
                                        double *krep, int *params) {
//    if(blocksnthreadse[1]>0) {
//

//        CUDAExclVolRepulsionenergy << < blocksnthreadse[0], blocksnthreadse[1],
//                (12 * blocksnthreadse[1]) * sizeof(double), stream >> >(coord, f, beadSet, krep, params, gU_i,
//                CUDAcommon::getCUDAvars().gculpritID,
//                CUDAcommon::getCUDAvars().gculpritFF,
//                CUDAcommon::getCUDAvars().gculpritinteraction, gFF, ginteraction);
////        CUDAcommon::handleerror(cudaEventRecord(event, stream));


//        CUDAcommon::handleerror(cudaGetLastError(),"CUDAExclVolRepulsionenergy", "CylinderExclVolumeRepulsion.cu");

//        double* gpu_Utot = CUDAcommon::getCUDAvars().gpu_energy;
//        addvector<<<1,1,0, stream>>>(gU_i,params, gU_sum, gpu_Utot);

//        auto cvars = CUDAcommon::getCUDAvars();
//        cvars.streamvec.push_back(&stream);
//        CUDAcommon::cudavars = cvars;

//        CUDAcommon::handleerror( cudaGetLastError() ,"CUDAExclVolRepulsionenergy", "CylinderExclVolumeRepulsion.cu");

//
//        return gU_sum;
//    }
//    else return NULL;
}

double* CylinderExclVolRepulsion::energy(double *coord, double *f, int *beadSet, double *krep, double *z, int *params) {

//    if(blocksnthreadse[1]>0) {
//        CUDAExclVolRepulsionenergy << < blocksnthreadse[0], blocksnthreadse[1],
//                (12 * blocksnthreadse[1]) * sizeof(double), stream >> >(coord, f, beadSet, krep, params, gU_i, z,
//                CUDAcommon::getCUDAvars().gculpritID,
//                CUDAcommon::getCUDAvars().gculpritFF,
//                CUDAcommon::getCUDAvars().gculpritinteraction, gFF, ginteraction);
//        CUDAcommon::handleerror(cudaGetLastError(),"CUDAExclVolRepulsionenergy", "CylinderExclVolumeRepulsion.cu");
//    }
    if(blocksnthreadsez[1]>0) {
        auto boolvarvec = CUDAcommon::cudavars.backtrackbools;
        CUDAExclVolRepulsionenergyz << < blocksnthreadsez[0], blocksnthreadsez[1],
                12 * blocksnthreadsez[1] * sizeof(double),stream >> > (coord, f, beadSet,
                krep, params, gU_i, CUDAcommon::cudavars.gpu_energyvec, z,
                CUDAcommon::getCUDAvars().gculpritID,
                CUDAcommon::getCUDAvars().gculpritFF,
                CUDAcommon::getCUDAvars().gculpritinteraction, gFF, ginteraction, boolvarvec.at(0),
                boolvarvec.at(1));
    }
    if(blocksnthreadse[1]<=0 && blocksnthreadsez[1]<=0)
        return NULL;
    else{
#ifdef CUDA_INDIVIDUAL_ESUM
        auto cvars = CUDAcommon::getCUDAvars();
        cvars.streamvec.push_back(&stream);
        CUDAcommon::cudavars = cvars;
        double* gpu_Utot = CUDAcommon::getCUDAvars().gpu_energy;
        resetdoublevariableCUDA<<<1,1,0,stream>>>(gU_sum);
        addvectorred2<<<bntaddvec2.at(2),bntaddvec2.at(3), bntaddvec2.at(3) * sizeof(double),stream>>>(gU_i,
                params, gU_sum, gpu_Utot);
#endif
                CUDAcommon::handleerror( cudaGetLastError() ,"CUDAExclVolRepulsionenergy",
"CylinderExclVolumeRepulsion.cu");
        return gU_sum;
    }
}

void CylinderExclVolRepulsion::forces(double *coord, double *f, int *beadSet, double *krep, int *params) {
//    cudaEvent_t start, stop;
//    CUDAcommon::handleerror(cudaEventCreate( &start));
//    CUDAcommon::handleerror(cudaEventCreate( &stop));
//    CUDAcommon::handleerror(cudaEventRecord( start, 0));

    if(blocksnthreadsf[1]>0) {
//        double *gU_ii;
//        double *gf1, *gf2, *gf3, *gf4, *gf5;
//        double *gc1, *gc2, *gcheckU;
//        double U_ii[blocksnthreadsf[0] * blocksnthreadsf[1]];
//        double c1[3 * blocksnthreadsf[0] * blocksnthreadsf[1]], c2[3 * blocksnthreadsf[0] * blocksnthreadsf[1]];
//        double F_i[3 * blocksnthreadsf[0] * blocksnthreadsf[1]];
//        double checkU[blocksnthreadsf[1]];

//        std::cout << "CEVF Number of Blocks: " << blocksnthreadsf[0] << endl;
//        std::cout << "Threads per block: " << blocksnthreadsf[1] << endl;

        //TODO  since the number of threads needed is constant through out the minimization, consider storing the pointer.
//        CUDAcommon::handleerror(cudaMalloc((void **) &gf1, 3 * blocksnthreadsf[0] * blocksnthreadsf[1] * sizeof
// (double)));
//        CUDAcommon::handleerror(cudaMalloc((void **) &gf2, 3 * blocksnthreadsf[0] * blocksnthreadsf[1] * sizeof
// (double)));
//        CUDAcommon::handleerror(cudaMalloc((void **) &gf3, 3 * blocksnthreadsf[0] * blocksnthreadsf[1] * sizeof
// (double)));
//        CUDAcommon::handleerror(cudaMalloc((void **) &gf4, 3 * blocksnthreadsf[0] * blocksnthreadsf[1] * sizeof
// (double)));
//        CUDAcommon::handleerror(
//                cudaMalloc((void **) &gf5, 45 * blocksnthreadsf[0] * blocksnthreadsf[1] * sizeof(double)));

//            size_t freeMem, totalMem;
//
//    cudaMemGetInfo(&freeMem, &totalMem);
//
//    std::cout<<"Memory "<<freeMem<<" "<<totalMem<<endl;
//        struct cudaDeviceProp properties;
//        cudaGetDeviceProperties(&properties, 0);
//        cout << "using " << properties.multiProcessorCount << " multiprocessors" << endl;
//        cout << "max threads per processor: " << properties.maxThreadsPerMultiProcessor << endl;
//        std::cout << 12 *blocksnthreadsf[1] * sizeof(double) << endl;
//        int blockSize;   // The launch configurator returned block size
//        int minGridSize; // The minimum grid size needed to achieve the
//        // maximum occupancy for a full device launch
//        int gridSize;    // The actual grid size needed, based on input size
//        cudaOccupancyMaxPotentialBlockSize( &minGridSize, &blockSize,
//                                            CUDAExclVolRepulsionforce, 0, 0);
//        gridSize = (blocksnthreadsf[0] * blocksnthreadsf[1] + blockSize -1) / blockSize;
//
//        std::cout<<gridSize<<" "<<blockSize<<endl;
//        size_t my_kernel_sm_size, count;

//        cudaOccupancyMaxPotentialBlockSizeVariableSMem( &minGridSize, &blockSize, CUDAExclVolRepulsionforce,
//                                                        my_kernel_sm_size, count);
//        std::cout<<minGridSize<<" "<<blockSize<<" "<<my_kernel_sm_size<<" "<<count<<endl;
//        int numblocks;
//        unsigned int flag;
//        cudaOccupancyMaxActiveBlocksPerMultiprocessorWithFlags(&numblocks,  CUDAExclVolRepulsionforce,blockSize,
//                                                               my_kernel_sm_size,flag);
//        std::cout<<numblocks<<" "<<blockSize<<" "<<my_kernel_sm_size<<endl;

        CUDAExclVolRepulsionforce << < blocksnthreadsf[0],blocksnthreadsf[1],
                12 *blocksnthreadsf[1] * sizeof(double),stream >> >
                (coord, f, beadSet, krep, params);

        CUDAcommon::handleerror(cudaGetLastError(),"CUDAExclVolRepulsionforce", "CylinderExclVolumeRepulsion.cu");
        auto cvars = CUDAcommon::getCUDAvars();
        cvars.streamvec.push_back(&stream);
        CUDAcommon::cudavars = cvars;
//    CUDAcommon::handleerror( cudaPeekAtLastError() );
       //CUDAcommon::handleerror(cudaDeviceSynchronize());

//    double f1 [3 * blocksnthreadsf[0]*blocksnthreadsf[1]];
//    double f2 [3 * blocksnthreadsf[0]*blocksnthreadsf[1]];
//    double f3 [3 * blocksnthreadsf[0]*blocksnthreadsf[1]];
//    double f4 [3 * blocksnthreadsf[0]*blocksnthreadsf[1]];
//    double f5 [45 * blocksnthreadsf[0]*blocksnthreadsf[1]];
//    cudaMemcpy(f1, gf1, 3 * blocksnthreadsf[0] * blocksnthreadsf[1] * sizeof(double), cudaMemcpyDeviceToHost);
//    cudaMemcpy(f2, gf2, 3 * blocksnthreadsf[0] * blocksnthreadsf[1] * sizeof(double), cudaMemcpyDeviceToHost);
//    cudaMemcpy(f3, gf3, 3 * blocksnthreadsf[0] * blocksnthreadsf[1] * sizeof(double), cudaMemcpyDeviceToHost);
//    cudaMemcpy(f4, gf4, 3 * blocksnthreadsf[0] * blocksnthreadsf[1] * sizeof(double), cudaMemcpyDeviceToHost);
//    cudaMemcpy(f5, gf5, 45 * blocksnthreadsf[0] * blocksnthreadsf[1] * sizeof(double), cudaMemcpyDeviceToHost);
//
//    for(auto i=0;i < blocksnthreadsf[0] * blocksnthreadsf[1]; i++) {
//        for(auto j=0;j<44;j++)
//            std::cout<<f5[44*i+j]<<" ";
//        std::cout << f1[3 * i] << " " << f1[3 * i + 1] << " " << f1[3 * i + 2] << " " << f2[3 * i] << " "
//                  << f2[3 * i + 1] << " " << f2[3 * i + 2] << " "
//                  << f3[3 * i] << " " << f3[3 * i + 1] << " " << f3[3 * i + 2] << " " << f4[3 * i] << " "
//                  << f4[3 * i + 1] << " " << f4[3 * i + 2] << endl;
//    }

//        CUDAcommon::handleerror(cudaFree(gf1));
//        CUDAcommon::handleerror(cudaFree(gf2));
//        CUDAcommon::handleerror(cudaFree(gf3));
//        CUDAcommon::handleerror(cudaFree(gf4));
//        CUDAcommon::handleerror(cudaFree(gf5));
    }
//    CUDAcommon::handleerror(cudaEventRecord( stop, 0));
//    CUDAcommon::handleerror(cudaEventSynchronize(stop));
//    float elaspedtime;
//    CUDAcommon::handleerror(cudaEventElapsedTime(&elaspedtime, start, stop));
//    CUDAvars cvars=CUDAcommon::getCUDAvars();
//    cvars.Ccforce += elaspedtime;
//    std::cout<<"C CFE "<<elaspedtime<<endl;
//    CUDAcommon::cudavars=cvars;
//    CUDAcommon::handleerror(cudaEventDestroy(start));
//    CUDAcommon::handleerror(cudaEventDestroy(stop));
}

void CylinderExclVolRepulsion::checkforculprit() {
    CUDAcommon::printculprit("Excluded Volume", "Cylinder Excluded Volume");
    int i = 0;
    cout<<"Printing culprit cylinders.."<<endl;
    for (auto cyl: Cylinder::getCylinders()) {
            auto id1 = cyl->getFirstBead()->getIndex();
            auto id2 = cyl->getSecondBead()->getIndex();
            if(id1 == CUDAcommon::getCUDAvars().culpritID[0] && id2 == CUDAcommon::getCUDAvars().culpritID[1])
                cyl->printSelf();
            else if(id1 == CUDAcommon::getCUDAvars().culpritID[2] && id2 == CUDAcommon::getCUDAvars().culpritID[3])
                cyl->printSelf();
        }
    exit(EXIT_FAILURE);
}

#endif
double CylinderExclVolRepulsion::energy(double *coord, int *beadSet, double *krep) {
    double *c1, *c2, *c3, *c4, *newc2, dsq, d, invDSquare, U, U_i;
    double a, b, c, e, F, AA, BB, CC, DD, EE, FF, GG, HH, JJ;
    double ATG1, ATG2, ATG3, ATG4;

    int nint = CylinderExclVolume<CylinderExclVolRepulsion>::numInteractions;
    int n = CylinderExclVolume<CylinderExclVolRepulsion>::n;

    U_i = 0.0;
    U = 0.0;
    newc2 = new double[3];
//    std::cout<<"SERL ecvol nint "<<nint<<endl;
    for (int i = 0; i < nint; i++) {

        c1 = &coord[3 * beadSet[n * i]];
        c2 = &coord[3 * beadSet[n * i + 1]];
        c3 = &coord[3 * beadSet[n * i + 2]];
        c4 = &coord[3 * beadSet[n * i + 3]];

//        if(areParallel(c1, c2, c3, c4)) std::cout<<"20"<<endl;
//        else if(areInPlane(c1, c2, c3, c4)) std::cout<<"11"<<endl;
//        else    std::cout<<"35007.0"<<endl;

        //check if parallel
        if(areParallel(c1, c2, c3, c4)) {
//            SysParams::exvolcounter[0] += 1;
            d = twoPointDistance(c1, c3);
            invDSquare =  1 / (d * d);
            U_i = krep[i] * invDSquare * invDSquare;
//            std::cout<<i<<" "<<U_i<<endl;
            if(fabs(U_i) == numeric_limits<double>::infinity()
               || U_i != U_i || U_i < -1.0) {

                //set culprit and return TODO
                return -1.0;
            }
            U += U_i;
            continue;
        }

        //check if in same plane
        if(areInPlane(c1, c2, c3, c4)) {
//            SysParams::exvolcounter[1] += 1;
            //slightly move point
            movePointOutOfPlane(c1, c2, c3, c4, newc2, 2, 0.01);
            c2 = newc2;
        }
//        else
//            SysParams::exvolcounter[2] += 1;


        a = scalarProduct(c1, c2, c1, c2);
        b = scalarProduct(c3, c4, c3, c4);
        c = scalarProduct(c3, c1, c3, c1);
        d = scalarProduct(c1, c2, c3, c4);
        e = scalarProduct(c1, c2, c3, c1);
        F = scalarProduct(c3, c4, c3, c1);

        AA = sqrt(a*c - e*e);
        BB = sqrt(b*c - F*F);

        CC = d*e - a*F;
        DD = b*e - d*F;

        EE = sqrt( a*(b + c - 2*F) - (d - e)*(d - e) );
        FF = sqrt( b*(a + c + 2*e) - (d + F)*(d + F) );

        GG = d*d - a*b - CC;
        HH = CC + GG - DD;
        JJ = c*(GG + CC) + e*DD - F*CC;


        ATG1 = atan( (a + e)/AA) - atan(e/AA);
        ATG2 = atan((a + e - d)/EE) - atan((e - d)/EE);
        ATG3 = atan((F)/BB) - atan((F - b)/BB);
        ATG4 = atan((d + F)/FF) - atan((d + F - b)/FF);

        U_i = 0.5 * krep[i]/ JJ * ( CC/AA*ATG1 + GG/EE*ATG2 + DD/BB*ATG3 + HH/FF*ATG4);
//        std::cout<<i<<" "<<U_i<<endl;

        if(fabs(U_i) == numeric_limits<double>::infinity()

           || U_i != U_i || U_i < -1.0) {

            //set culprit and return TODO

            return -1;
        }
        U += U_i;
    }
    delete [] newc2;
//    std::cout<<"total energy serial "<<U<<endl;
    return U;

}


double CylinderExclVolRepulsion::energy(double *coord, double *f, int *beadSet,
                                        double *krep, double z) {


    double d, dsq, invDSquare, U, U_i, *f1, *f2, *f3, *f4;
    double a, b, c, e, F, AA, BB, CC, DD, EE, FF, GG, HH, JJ;
    double ATG1, ATG2, ATG3, ATG4;
    double *c1us, *c2us, *c3us, *c4us;
    double *c1 = new double[3];//stretched
    double *c2 = new double[3];
    double *c3 = new double[3];
    double *c4 = new double[3];
//    double *c1us = new double[3];//unstretched
//    double *c2us = new double[3];
//    double *c3us = new double[3];
//    double *c4us = new double[3];
    double *newc2 = new double[3];

    int nint = CylinderExclVolume<CylinderExclVolRepulsion>::numInteractions;
    int n = CylinderExclVolume<CylinderExclVolRepulsion>::n;

    U_i = 0.0;
    U = 0.0;
//std::cout<<"-----------"<<endl;
    for (int i = 0; i < nint; i++) {

        c1us = &coord[3 * beadSet[n * i]];
        c2us = &coord[3 * beadSet[n * i +1]];
        c3us = &coord[3 * beadSet[n * i +2]];
        c4us = &coord[3 * beadSet[n * i +3]];

//        memcpy(c1us, &coord[3 * beadSet[n * i]], 3 * sizeof(double));
//        memcpy(c2us, &coord[3 * beadSet[n * i + 1]], 3 * sizeof(double));
//        memcpy(c3us, &coord[3 * beadSet[n * i + 2]], 3 * sizeof(double));
//        memcpy(c4us, &coord[3 * beadSet[n * i + 3]], 3 * sizeof(double));

        //stretch coords
        f1 = &f[3 * beadSet[n * i]];
        f2 = &f[3 * beadSet[n * i + 1]];
        f3 = &f[3 * beadSet[n * i + 2]];
        f4 = &f[3 * beadSet[n * i + 3]];

//        cout.precision(dbl::max_digits10);
//        std::cout<<i<<" BEFORE "<<c1us[0]<<" "<<c1us[1]<<" "<<c1us[2]<<" "<<c2us[0]<<" "<<c2us[1]<<" "
//                ""<<c2us[2]<<" "<<c3us[0]<<" "
//                ""<<c3us[1]<<" "
//                ""<<c3us[2]<<" "<<c4us[0]<<" "<<c4us[1]<<" "<<c4us[2]<<endl;
//        std::cout<<"Force "<<f1[0]<<" "
//                ""<<f1[1]<<" "<<f1[2]<<" "<<f2[0]<<" "
//                ""<<f2[1]<<" "<<f2[2]<<" "<<f3[0]<<" "<<f3[1]<<" "<<f3[2]<<" "<<f4[0]<<" "<<f4[1]<<" "<<f4[2]<<endl;
        c1[0] = c1us[0] + z * f1[0];
        c1[1] = c1us[1] + z * f1[1];
        c1[2] = c1us[2] + z * f1[2];
        c2[0] = c2us[0] + z * f2[0];
        c2[1] = c2us[1] + z * f2[1];
        c2[2] = c2us[2] + z * f2[2];
        c3[0] = c3us[0] + z * f3[0];
        c3[1] = c3us[1] + z * f3[1];
        c3[2] = c3us[2] + z * f3[2];
        c4[0] = c4us[0] + z * f4[0];
        c4[1] = c4us[1] + z * f4[1];
        c4[2] = c4us[2] + z * f4[2];
//        std::cout<<"AFTER "<<c1[0]<<" "<<c1[1]<<" "<<c1[2]<<" "<<c2[0]<<" "<<c2[1]<<" "
//                ""<<c2[1]<<" "<<
//                 c3[0]<<" "<<c3[1]<<" "<<c3[2]<<" "<<c4[0]<<" "<<c4[1]<<" "<<c4[2]<<endl;
        //check if parallel
        if(areParallel(c1, c2, c3, c4)) {
//            SysParams::exvolcounterz[0] += 1;
            d = twoPointDistance(c1, c3);
            invDSquare =  1 / (d * d);
            U_i = krep[i] * invDSquare * invDSquare;
//            std::cout<<"P Energy"<<U_i<<endl;
            if(fabs(U_i) == numeric_limits<double>::infinity()
               || U_i != U_i || U_i < -1.0) {

                //set culprit and return TODO
                return -1;
            }
            U += U_i;
            continue;
        }

        //check if in same plane
        if(areInPlane(c1, c2, c3, c4)) {
//            SysParams::exvolcounterz[1] += 1;
            //slightly move point
            movePointOutOfPlane(c1, c2, c3, c4, newc2, 2, 0.01);
            c2 = newc2;
            std::cout<<"move"<<endl;
//            std::cout<<i<<" 2.0 "<<c1[0]<<" "<<c1[1]<<" "<<c1[2]<<" "<<c2[0]<<" "<<c2[1]<<" "<<c2[2]<<" "<<c3[0]<<" "
//                    ""<<c3[1]<<" "
//                             ""<<c3[2]<<" "<<c4[0]<<" "<<c4[1]<<" "<<c4[2]<<" "<<U_i<<endl;
        }
//        else{
//            SysParams::exvolcounterz[2] += 1;
//        }
        a = scalarProduct(c1, c2, c1, c2);
        b = scalarProduct(c3, c4, c3, c4);
        c = scalarProduct(c3, c1, c3, c1);
        d = scalarProduct(c1, c2, c3, c4);
        e = scalarProduct(c1, c2, c3, c1);
        F = scalarProduct(c3, c4, c3, c1);
//    std::cout<<i<<" "<<a<<" "<<b<<" "<<c<<" "<<d<<" "<<e<<" "<<F<<endl;
        AA = sqrt(a*c - e*e);
        BB = sqrt(b*c - F*F);

        CC = d*e - a*F;
        DD = b*e - d*F;

        EE = sqrt( a*(b + c - 2*F) - (d - e)*(d - e) );
        FF = sqrt( b*(a + c + 2*e) - (d + F)*(d + F) );

        GG = d*d - a*b - CC;
        HH = CC + GG - DD;
        JJ = c*(GG + CC) + e*DD - F*CC;


        ATG1 = atan( (a + e)/AA) - atan(e/AA);
        ATG2 = atan((a + e - d)/EE) - atan((e - d)/EE);
        ATG3 = atan((F)/BB) - atan((F - b)/BB);
        ATG4 = atan((d + F)/FF) - atan((d + F - b)/FF);

        U_i = 0.5 * krep[i]/ JJ * ( CC/AA*ATG1 + GG/EE*ATG2 + DD/BB*ATG3 + HH/FF*ATG4);
//        std::cout<<"N energy "<<U_i<<endl;
        if(fabs(U_i) == numeric_limits<double>::infinity()
           || U_i != U_i || U_i < -1.0) {

            //set culprit and return TODO

            return -1.0;
        }
//        std::cout<<i<<" "<<U_i<<" "<<a<<" "<<b<<" "<<c<<" "<<d<<" "<<e<<" "<<F<<" "<<AA<<" "<<BB<<" "<<CC<<" "<<DD<<" "
//                ""<<EE<<" "
//                ""<<FF<<""
//                " "<<GG<<" "<<HH<<" "<<JJ<<" "<<ATG1<<" "<<ATG2<<" "<<ATG3<<" "<<ATG4<<endl;
//        std::cout<<U_i<<endl;
        U += U_i;
    }
//    std::cout<<"Total energy serial "<<U<<endl;
    delete [] c1;
    delete [] c2;
    delete [] c3;
    delete [] c4;
    delete [] newc2;
//    delete [] c1us;
//    delete [] c2us;
//    delete [] c3us;
//    delete [] c4us;
    return U;
}

void CylinderExclVolRepulsion::forces(double *coord, double *f, int *beadSet, double *krep) {
//cout.precision(10);
//    clock_t start, stop;
//    float elapsedtime;
//    start = clock();

//    cout.precision(dbl::max_digits10); //TODO remove precision.
    double *c1, *c2, *c3, *c4, *newc2, d, dsq, invDSquare, U, *f1, *f2, *f3, *f4;
    double a, b, c, e, F, AA, BB, CC, DD, EE, FF, GG, HH, JJ, invJJ;
    double ATG1, ATG2, ATG3, ATG4;
    double A1, A2, E1, E2, B1, B2, F1, F2, A11, A12, A13, A14;
    double E11, E12, E13, E14, B11, B12, B13, B14, F11, F12, F13, F14;

    int nint = CylinderExclVolume<CylinderExclVolRepulsion>::numInteractions;
    int n = CylinderExclVolume<CylinderExclVolRepulsion>::n;

//            for(auto i=0;i<CGMethod::N;i++)
//            std::cout<<f[i]<<" ";
//        std::cout<<endl;
    //std::cout<<"Excl vol nint "<<nint<<endl;
    newc2 = new double[3];
    for (int i = 0; i < nint; i++) {
//        std::cout<<beadSet[n * i]<<" "<<beadSet[n * i+1]<<" "<<beadSet[n * i+2]<<" "<<beadSet[n * i+3]<<endl;
        c1 = &coord[3 * beadSet[n * i]];
        c2 = &coord[3 * beadSet[n * i + 1]];
        c3 = &coord[3 * beadSet[n * i + 2]];
        c4 = &coord[3 * beadSet[n * i + 3]];

        //stretch coords
        f1 = &f[3 * beadSet[n * i]];
        f2 = &f[3 * beadSet[n * i + 1]];
        f3 = &f[3 * beadSet[n * i + 2]];
        f4 = &f[3 * beadSet[n * i + 3]];

//    std::cout<<c1[0]<<" "<<c1[1]<<" "<<c1[2]<<" "<<
//         c2[0]<<" "<<c2[1]<<" "<<c2[2]<<" "<<
//         c3[0]<<" "<<c3[1]<<" "<<c3[2]<<" "<<
//         c4[0]<<" "<<c4[1]<<" "<<c4[2]<<endl;
        //check if parallel
        if(areParallel(c1, c2, c3, c4)) {

            dsq = twoPointDistancesquared(c1, c3);
            invDSquare =  1 / dsq;
            U = krep[i] * invDSquare * invDSquare;

            double f0 = 4 * krep[i] * invDSquare * invDSquare * invDSquare;

            f1[0] += - f0 * (c3[0] - c1[0]);
            f1[1] += - f0 * (c3[1] - c1[1]);
            f1[2] += - f0 * (c3[2] - c1[2]);

            f2[0] += - f0 * (c4[0] - c2[0]);
            f2[1] += - f0 * (c4[1] - c2[1]);
            f2[2] += - f0 * (c4[2] - c2[2]);

            f3[0] += f0 * (c3[0] - c1[0]);
            f3[1] += f0 * (c3[1] - c1[1]);
            f3[2] += f0 * (c3[2] - c1[2]);

            f4[0] += f0 * (c4[0] - c2[0]);
            f4[1] += f0 * (c4[1] - c2[1]);
            f4[2] += f0 * (c4[2] - c2[2]);
#ifdef DETAILEDOUTPUT
            std::cout<<"P "<<f1[0]<<" "<<f1[1]<<" "<<f1[2]<<" "<<
                            f2[0]<<" "<<f2[1]<<" "<<f2[2]<<" "<<
                            f3[0]<<" "<<f3[1]<<" "<<f3[2]<<" "<<
                            f4[0]<<" "<<f4[1]<<" "<<f4[2]<<endl;
#endif
            continue;
        }

        //check if in same plane
        if(areInPlane(c1, c2, c3, c4)) {

            //slightly move point
            movePointOutOfPlane(c1, c2, c3, c4, newc2, 2, 0.01);
            c2 = newc2;
//            delete c2;
#ifdef DETAILEDOUTPUT
            std::cout<<"Mv"<<c1[0]<<" "<<c1[1]<<" "<<c1[2]<<" "<<
                     c2[0]<<" "<<c2[1]<<" "<<c2[2]<<" "<<
                     c3[0]<<" "<<c3[1]<<" "<<c3[2]<<" "<<
                     c4[0]<<" "<<c4[1]<<" "<<c4[2]<<endl;
            std::cout<<"M ";
#endif
        }
#ifdef DETAILEDOUTPUT
        else{
            std::cout<<"N ";
        }
#endif
        a = scalarProduct(c1, c2, c1, c2);
        b = scalarProduct(c3, c4, c3, c4);
        c = scalarProduct(c3, c1, c3, c1);
        d = scalarProduct(c1, c2, c3, c4);
        e = scalarProduct(c1, c2, c3, c1);
        F = scalarProduct(c3, c4, c3, c1);
//        std::cout<<"N "<<a<<" "<<b<<" "<<c<<" "<<d<<" "<<e<<" "<<F<<endl;
        AA = sqrt(a*c - e*e);
        BB = sqrt(b*c - F*F);

        CC = d*e - a*F;
        DD = b*e - d*F;

        EE = sqrt( a*(b + c - 2*F) - (d - e)*(d - e) );
        FF = sqrt( b*(a + c + 2*e) - (d + F)*(d + F) );

        GG = d*d - a*b - CC;
        HH = CC + GG - DD;
        JJ = c*(GG + CC) + e*DD - F*CC;
//        std::cout<<"N2 "<<AA<<" "<<BB<<" "<<CC<<" "<<DD<<" "<<EE<<" "<<FF<<" "<<GG<<" "<<HH<<" "<<JJ<<endl;
        invJJ = 1/JJ;

        ATG1 = atan( (a + e)/AA) - atan(e/AA);
        ATG2 = atan((a + e - d)/EE) - atan((e - d)/EE);
        ATG3 = atan((F)/BB) - atan((F - b)/BB);
        ATG4 = atan((d + F)/FF) - atan((d + F - b)/FF);
//        std::cout<<"N3 "<<ATG1<<" "<<ATG2<<" "<<ATG3<<" "<<ATG4<<endl;
//         U = 0.5 * krep[i]/ JJ * ( CC/AA*ATG1 + GG/EE*ATG2 + DD/BB*ATG3 + HH/FF*ATG4);
#ifdef DETAILEDOUTPUT
        std::cout<<a<<" "<<b<<" "<<c<<" "<<d<<" "<<e<<" "<<F<<" "<<AA<<" "<<BB<<" "<<CC<<" "
                ""<<DD<<" "<<EE<<" "<<FF<<" "<<GG<<" "<<HH<<" "<<JJ<<" "<<ATG1<<" "<<ATG2<<" "
                         ""<<ATG3<<" "<<ATG4<<" "<<U<<" "<<krep[i]<<endl;
#endif
        U = 0.5 * krep[i]*invJJ * ( CC/AA*ATG1 + GG/EE*ATG2 + DD/BB*ATG3 + HH/FF*ATG4);
//        std::cout<<U<<endl;
        A1 = AA*AA/(AA*AA + e*e);
        A2 = AA*AA/(AA*AA + (a + e)*(a + e));

        E1 = EE*EE/(EE*EE + (a + e - d)*(a + e - d));
        E2 = EE*EE/(EE*EE + (e - d)*(e - d));

        B1 = BB*BB/(BB*BB + (F - b)*(F - b));
        B2 = BB*BB/(BB*BB + F*F);

        F1 = FF*FF/(FF*FF + (d + F - b)*(d + F - b));
        F2 = FF*FF/(FF*FF + (d + F)*(d + F));
//        std::cout<<"N4 "<<U<<" "<<krep[i]<<" "<<A1<<" "<<A2<<" "<<E1<<" "<<E2<<" "<<B1<<" "<<B2<<" "<<F1<<" "<<F2<<endl;
        A11 = ATG1/AA;
        A12 = -((ATG1*CC)/(AA*AA)) + (A1*CC*e)/(AA*AA*AA) -
              (A2*CC*(a + e))/(AA*AA*AA);
        A13 = -((A1*CC)/(AA*AA)) + (A2*CC)/(AA*AA);
        A14 = (A2*CC)/(AA*AA);

        E11 = ATG2/EE;
        E12 = (E2*(-a + d - e)*GG)/(EE*EE*EE) + (E1*(-d + e)*GG)/(EE*EE*EE) -
              (ATG2*GG)/(EE*EE);
        E13 = -((E1*GG)/(EE*EE)) + (E2*GG)/(EE*EE);
        E14 = (E2*GG)/(EE*EE);

        B11 = ATG3/BB;
        B12 = -((ATG3*DD)/(BB*BB)) - (B2*DD*F)/(BB*BB*BB) +
              (B1*DD*(-b + F))/(BB*BB*BB);
        B13 = -((B1*DD)/(BB*BB)) + (B2*DD)/(BB*BB);
        B14 = (B1*DD)/(BB*BB);

        F11 = ATG4/FF;
        F12 = (F2*(-d - F)*HH)/(FF*FF*FF) + (F1*(-b + d + F)*HH)/(FF*FF*FF) -
              (ATG4*HH)/(FF*FF);
        F13 = -((F1*HH)/(FF*FF)) + (F2*HH)/(FF*FF);
        F14 = (F1*HH)/(FF*FF);
//        std::cout<<A1<<" "<<A2<<" "<<E1<<" "<<E2<<" "<<B1<<" "<<B2<<" "<<F1<<" "<<F2<<" "
//                ""<<A11<<" "<<A12<<" "<<A13<<" "<<A14<<" "<<E11<<" "<<E12<<" "<<E13<<" "
//                         ""<<E14<<" "<<B11<<" "<<B12<<" "<<B13<<" "<<B14<<" "<<F11<<" "<<F12<<" "
//                         ""<<F13<<" "<<F14<<endl;
        f1[0] +=  - 0.5*invJJ*( (c2[0] - c1[0] ) *( A13 + E13 + B11*b - F11*b + A11*d - E11*d - 2*U*b*e - (A12*e)/AA + (E12*(d - e))/EE + 2*U*d*F - 2*U*(b*e - d*F) + (F12*b)/FF - 2*(A14 + E14 - E11*b - F11*b + 2*U*b*c + (A12*c)/(2*AA) + (E12*(b + c - 2*F))/(2*EE) - A11*F + E11*F - 2*U*F*F + (F12*b)/(2*FF)) ) + (c4[0] - c3[0] ) *(B13 + E13 - A11*a + E11*a - B11*d - 2*E11*d - F11*d + 4*U*c*d - A11*e + E11*e + 2*U*d*e - (E12*a)/EE + (E12*(d - e))/EE + B11*F - F11*F - 2*U*a*F - 4*U*e*F + 2*U*(d*e - a*F) - (B12*F)/BB) +  (c1[0] - c3[0] )* (-A13 - E13 - B11*b + F11*b - A11*d + E11*d + 2*U*b*e + (A12*e)/AA - (E12*(d - e))/EE - 2*U*d*F + 2*U*(b*e - d*F) - (F12*b)/FF + 2*(-2*U*((-a)*b + d*d) + (A12*a)/(2*AA) + (E12*a)/(2*EE) +(B12*b)/(2*BB) + (F12*b)/(2*FF))) );

        f1[1] +=  - 0.5*invJJ*( (c2[1] - c1[1] ) *( A13 + E13 + B11*b - F11*b + A11*d - E11*d - 2*U*b*e - (A12*e)/AA + (E12*(d - e))/EE + 2*U*d*F - 2*U*(b*e - d*F) + (F12*b)/FF - 2*(A14 + E14 - E11*b - F11*b + 2*U*b*c + (A12*c)/(2*AA) + (E12*(b + c - 2*F))/(2*EE) - A11*F + E11*F - 2*U*F*F + (F12*b)/(2*FF)) ) + (c4[1] - c3[1] ) *(B13 + E13 - A11*a + E11*a - B11*d - 2*E11*d - F11*d + 4*U*c*d - A11*e + E11*e + 2*U*d*e - (E12*a)/EE + (E12*(d - e))/EE + B11*F - F11*F - 2*U*a*F - 4*U*e*F + 2*U*(d*e - a*F) - (B12*F)/BB) +  (c1[1] - c3[1] )* (-A13 - E13 - B11*b + F11*b - A11*d + E11*d + 2*U*b*e + (A12*e)/AA - (E12*(d - e))/EE - 2*U*d*F + 2*U*(b*e - d*F) - (F12*b)/FF + 2*(-2*U*((-a)*b + d*d) + (A12*a)/(2*AA) + (E12*a)/(2*EE) +(B12*b)/(2*BB) + (F12*b)/(2*FF))) );

        f1[2] +=  - 0.5*invJJ*( (c2[2] - c1[2] ) *( A13 + E13 + B11*b - F11*b + A11*d - E11*d - 2*U*b*e - (A12*e)/AA + (E12*(d - e))/EE + 2*U*d*F - 2*U*(b*e - d*F) + (F12*b)/FF - 2*(A14 + E14 - E11*b - F11*b + 2*U*b*c + (A12*c)/(2*AA) + (E12*(b + c - 2*F))/(2*EE) - A11*F + E11*F - 2*U*F*F + (F12*b)/(2*FF)) ) + (c4[2] - c3[2] ) *(B13 + E13 - A11*a + E11*a - B11*d - 2*E11*d - F11*d + 4*U*c*d - A11*e + E11*e + 2*U*d*e - (E12*a)/EE + (E12*(d - e))/EE + B11*F - F11*F - 2*U*a*F - 4*U*e*F + 2*U*(d*e - a*F) - (B12*F)/BB) +  (c1[2] - c3[2] )* (-A13 - E13 - B11*b + F11*b - A11*d + E11*d + 2*U*b*e + (A12*e)/AA - (E12*(d - e))/EE - 2*U*d*F + 2*U*(b*e - d*F) - (F12*b)/FF + 2*(-2*U*((-a)*b + d*d) + (A12*a)/(2*AA) + (E12*a)/(2*EE) +(B12*b)/(2*BB) + (F12*b)/(2*FF))) );


        f2[0] +=  - invJJ*( (c2[0] - c1[0] )*( A14+E14-E11*b-F11*b+2*U*b*c+(A12*c)/(2*AA)+(E12*(b+c-2*F))/(2*EE)-A11*F+E11*F-2*U*F*F+(F12*b)/(2*FF) ) + 0.5*(c4[0] - c3[0])*(-E13 + F13 + 2*E11*d + 2*F11*d - 4*U*c*d + A11*e - E11*e - (E12*(d - e))/EE - B11*F + F11*F + 4*U*e*F - (F12*(d + F))/FF)  + 0.5*(c1[0] - c3[0] )* (A13 + E13 + B11*b - F11*b + A11*d - E11*d - 2*U*b*e - (A12*e)/AA + (E12*(d - e))/EE + 2*U*d*F - 2*U*(b*e - d*F) + (F12*b)/FF) );

        f2[1] += - invJJ*( (c2[1] - c1[1] )*( A14+E14-E11*b-F11*b+2*U*b*c+(A12*c)/(2*AA)+(E12*(b+c-2*F))/(2*EE)-A11*F+E11*F-2*U*F*F+(F12*b)/(2*FF) ) + 0.5*(c4[1] - c3[1])*(-E13 + F13 + 2*E11*d + 2*F11*d - 4*U*c*d + A11*e - E11*e - (E12*(d - e))/EE - B11*F + F11*F + 4*U*e*F - (F12*(d + F))/FF)  + 0.5*(c1[1] - c3[1] )* (A13 + E13 + B11*b - F11*b + A11*d - E11*d - 2*U*b*e - (A12*e)/AA + (E12*(d - e))/EE + 2*U*d*F - 2*U*(b*e - d*F) + (F12*b)/FF) );

        f2[2] += - invJJ*( (c2[2] - c1[2] )*( A14+E14-E11*b-F11*b+2*U*b*c+(A12*c)/(2*AA)+(E12*(b+c-2*F))/(2*EE)-A11*F+E11*F-2*U*F*F+(F12*b)/(2*FF) ) + 0.5*(c4[2] - c3[2])*(-E13 + F13 + 2*E11*d + 2*F11*d - 4*U*c*d + A11*e - E11*e - (E12*(d - e))/EE - B11*F + F11*F + 4*U*e*F - (F12*(d + F))/FF)  + 0.5*(c1[2] - c3[2] )* (A13 + E13 + B11*b - F11*b + A11*d - E11*d - 2*U*b*e - (A12*e)/AA + (E12*(d - e))/EE + 2*U*d*F - 2*U*(b*e - d*F) + (F12*b)/FF) );

        f3[0] +=  - 0.5*invJJ*( (c2[0] - c1[0] )*(-A13 - F13 - B11*b + F11*b - A11*d - E11*d - 2*F11*d + 4*U*c*d - A11*e + E11*e + 2*U*b*e + (A12*e)/AA + B11*F - F11*F - 2*U*d*F - 4*U*e*F + 2*U*(b*e - d*F) - (F12*b)/FF + (F12*(d + F))/FF) + (c4[0] - c3[0] )*(-B13 - F13 + A11*a - E11*a + B11*d - F11*d - 2*U*d*e + (E12*a)/EE + 2*U*a*F - 2*U*(d*e - a*F) + (B12*F)/BB + (F12*(d + F))/FF - 2*(B14 + F14 - E11*a - F11*a + 2*U*a*c + B11*e - F11*e - 2*U*e*e + (E12*a)/(2*EE) + (B12*c)/(2*BB) + (F12*(a + c + 2*e))/(2*FF))) + (c1[0] - c3[0] ) * (-B13 - F13 + A11*a - E11*a + B11*d - F11*d - 2*U*d*e + (E12*a)/EE + 2*U*a*F - 2*U*(d*e - a*F) + (B12*F)/BB + (F12*(d + F))/FF - 2*(-2*U*((-a)*b + d*d) + (A12*a)/(2*AA) + (E12* a)/(2*EE) + (B12*b)/(2*BB) + (F12*b)/(2*FF))) );

        f3[1] +=  - 0.5*invJJ*( (c2[1] - c1[1] )*(-A13 - F13 - B11*b + F11*b - A11*d - E11*d - 2*F11*d + 4*U*c*d - A11*e + E11*e + 2*U*b*e + (A12*e)/AA + B11*F - F11*F - 2*U*d*F - 4*U*e*F + 2*U*(b*e - d*F) - (F12*b)/FF + (F12*(d + F))/FF) + (c4[1] - c3[1] )*(-B13 - F13 + A11*a - E11*a + B11*d - F11*d - 2*U*d*e + (E12*a)/EE + 2*U*a*F - 2*U*(d*e - a*F) + (B12*F)/BB + (F12*(d + F))/FF - 2*(B14 + F14 - E11*a - F11*a + 2*U*a*c + B11*e - F11*e - 2*U*e*e + (E12*a)/(2*EE) + (B12*c)/(2*BB) + (F12*(a + c + 2*e))/(2*FF))) + (c1[1] - c3[1] ) * (-B13 - F13 + A11*a - E11*a + B11*d - F11*d - 2*U*d*e + (E12*a)/EE + 2*U*a*F - 2*U*(d*e - a*F) + (B12*F)/BB + (F12*(d + F))/FF - 2*(-2*U*((-a)*b + d*d) + (A12*a)/(2*AA) + (E12* a)/(2*EE) + (B12*b)/(2*BB) + (F12*b)/(2*FF))) ) ;

        f3[2] +=  - 0.5*invJJ*( (c2[2] - c1[2] )*(-A13 - F13 - B11*b + F11*b - A11*d - E11*d - 2*F11*d + 4*U*c*d - A11*e + E11*e + 2*U*b*e + (A12*e)/AA + B11*F - F11*F - 2*U*d*F - 4*U*e*F + 2*U*(b*e - d*F) - (F12*b)/FF + (F12*(d + F))/FF) + (c4[2] - c3[2] )*(-B13 - F13 + A11*a - E11*a + B11*d - F11*d - 2*U*d*e + (E12*a)/EE + 2*U*a*F - 2*U*(d*e - a*F) + (B12*F)/BB + (F12*(d + F))/FF - 2*(B14 + F14 - E11*a - F11*a + 2*U*a*c + B11*e - F11*e - 2*U*e*e + (E12*a)/(2*EE) + (B12*c)/(2*BB) + (F12*(a + c + 2*e))/(2*FF))) + (c1[2] - c3[2] ) * (-B13 - F13 + A11*a - E11*a + B11*d - F11*d - 2*U*d*e + (E12*a)/EE + 2*U*a*F - 2*U*(d*e - a*F) + (B12*F)/BB + (F12*(d + F))/FF - 2*(-2*U*((-a)*b + d*d) + (A12*a)/(2*AA) + (E12* a)/(2*EE) + (B12*b)/(2*BB) + (F12*b)/(2*FF))) );


        f4[0] +=  - invJJ*( 0.5*(c2[0] - c1[0] )*( -E13 + F13 + 2*E11*d + 2*F11*d - 4*U*c*d + A11*e - E11*e - (E12*(d - e))/EE - B11*F + F11*F +4*U*e*F - (F12*(d + F))/FF ) + (c4[0] - c3[0])*(B14 + F14 - E11*a - F11*a + 2*U*a*c + B11*e - F11*e - 2*U*e*e + (E12*a)/(2*EE) + (B12*c)/(2*BB) + (F12*(a + c + 2*e))/(2*FF))  + 0.5*(c1[0] - c3[0] )* (B13 + F13 - A11*a + E11*a - B11*d + F11*d + 2*U*d*e - (E12*a)/EE - 2*U*a*F + 2*U*(d*e - a*F) - (B12*F)/BB - (F12*(d + F))/FF) )  ;

        f4[1] +=  - invJJ*( 0.5*(c2[1] - c1[1] )*( -E13 + F13 + 2*E11*d + 2*F11*d - 4*U*c*d + A11*e - E11*e - (E12*(d - e))/EE - B11*F + F11*F +4*U*e*F - (F12*(d + F))/FF ) + (c4[1] - c3[1])*(B14 + F14 - E11*a - F11*a + 2*U*a*c + B11*e - F11*e - 2*U*e*e + (E12*a)/(2*EE) + (B12*c)/(2*BB) + (F12*(a + c + 2*e))/(2*FF))  + 0.5*(c1[1] - c3[1] )* (B13 + F13 - A11*a + E11*a - B11*d + F11*d + 2*U*d*e - (E12*a)/EE - 2*U*a*F + 2*U*(d*e - a*F) - (B12*F)/BB - (F12*(d + F))/FF) ) ;

        f4[2] +=  - invJJ*( 0.5*(c2[2] - c1[2] )*( -E13 + F13 + 2*E11*d + 2*F11*d - 4*U*c*d + A11*e - E11*e - (E12*(d - e))/EE - B11*F + F11*F +4*U*e*F - (F12*(d + F))/FF ) + (c4[2] - c3[2])*(B14 + F14 - E11*a - F11*a + 2*U*a*c + B11*e - F11*e - 2*U*e*e + (E12*a)/(2*EE) + (B12*c)/(2*BB) + (F12*(a + c + 2*e))/(2*FF))  + 0.5*(c1[2] - c3[2] )* (B13 + F13 - A11*a + E11*a - B11*d + F11*d + 2*U*d*e - (E12*a)/EE - 2*U*a*F + 2*U*(d*e - a*F) - (B12*F)/BB - (F12*(d + F))/FF) ) ;


//        double fc1[3], fc2[3], fc3[3], fc4[3];
//        fc1[0] =  - 0.5*invJJ*( (c2[0] - c1[0] ) *( A13 + E13 + B11*b - F11*b + A11*d - E11*d - 2*U*b*e - (A12*e)/AA + (E12*(d - e))/EE + 2*U*d*F - 2*U*(b*e - d*F) + (F12*b)/FF - 2*(A14 + E14 - E11*b - F11*b + 2*U*b*c + (A12*c)/(2*AA) + (E12*(b + c - 2*F))/(2*EE) - A11*F + E11*F - 2*U*F*F + (F12*b)/(2*FF)) ) + (c4[0] - c3[0] ) *(B13 + E13 - A11*a + E11*a - B11*d - 2*E11*d - F11*d + 4*U*c*d - A11*e + E11*e + 2*U*d*e - (E12*a)/EE + (E12*(d - e))/EE + B11*F - F11*F - 2*U*a*F - 4*U*e*F + 2*U*(d*e - a*F) - (B12*F)/BB) +  (c1[0] - c3[0] )* (-A13 - E13 - B11*b + F11*b - A11*d + E11*d + 2*U*b*e + (A12*e)/AA - (E12*(d - e))/EE - 2*U*d*F + 2*U*(b*e - d*F) - (F12*b)/FF + 2*(-2*U*((-a)*b + d*d) + (A12*a)/(2*AA) + (E12*a)/(2*EE) +(B12*b)/(2*BB) + (F12*b)/(2*FF))) );
//
//        fc1[1] =  - 0.5*invJJ*( (c2[1] - c1[1] ) *( A13 + E13 + B11*b - F11*b + A11*d - E11*d - 2*U*b*e - (A12*e)/AA + (E12*(d - e))/EE + 2*U*d*F - 2*U*(b*e - d*F) + (F12*b)/FF - 2*(A14 + E14 - E11*b - F11*b + 2*U*b*c + (A12*c)/(2*AA) + (E12*(b + c - 2*F))/(2*EE) - A11*F + E11*F - 2*U*F*F + (F12*b)/(2*FF)) ) + (c4[1] - c3[1] ) *(B13 + E13 - A11*a + E11*a - B11*d - 2*E11*d - F11*d + 4*U*c*d - A11*e + E11*e + 2*U*d*e - (E12*a)/EE + (E12*(d - e))/EE + B11*F - F11*F - 2*U*a*F - 4*U*e*F + 2*U*(d*e - a*F) - (B12*F)/BB) +  (c1[1] - c3[1] )* (-A13 - E13 - B11*b + F11*b - A11*d + E11*d + 2*U*b*e + (A12*e)/AA - (E12*(d - e))/EE - 2*U*d*F + 2*U*(b*e - d*F) - (F12*b)/FF + 2*(-2*U*((-a)*b + d*d) + (A12*a)/(2*AA) + (E12*a)/(2*EE) +(B12*b)/(2*BB) + (F12*b)/(2*FF))) );
//
//        fc1[2] =  - 0.5*invJJ*( (c2[2] - c1[2] ) *( A13 + E13 + B11*b - F11*b + A11*d - E11*d - 2*U*b*e - (A12*e)/AA + (E12*(d - e))/EE + 2*U*d*F - 2*U*(b*e - d*F) + (F12*b)/FF - 2*(A14 + E14 - E11*b - F11*b + 2*U*b*c + (A12*c)/(2*AA) + (E12*(b + c - 2*F))/(2*EE) - A11*F + E11*F - 2*U*F*F + (F12*b)/(2*FF)) ) + (c4[2] - c3[2] ) *(B13 + E13 - A11*a + E11*a - B11*d - 2*E11*d - F11*d + 4*U*c*d - A11*e + E11*e + 2*U*d*e - (E12*a)/EE + (E12*(d - e))/EE + B11*F - F11*F - 2*U*a*F - 4*U*e*F + 2*U*(d*e - a*F) - (B12*F)/BB) +  (c1[2] - c3[2] )* (-A13 - E13 - B11*b + F11*b - A11*d + E11*d + 2*U*b*e + (A12*e)/AA - (E12*(d - e))/EE - 2*U*d*F + 2*U*(b*e - d*F) - (F12*b)/FF + 2*(-2*U*((-a)*b + d*d) + (A12*a)/(2*AA) + (E12*a)/(2*EE) +(B12*b)/(2*BB) + (F12*b)/(2*FF))) );
//
//
//        fc2[0] =  - invJJ*( (c2[0] - c1[0] )*( A14+E14-E11*b-F11*b+2*U*b*c+(A12*c)/(2*AA)+(E12*(b+c-2*F))/(2*EE)-A11*F+E11*F-2*U*F*F+(F12*b)/(2*FF) ) + 0.5*(c4[0] - c3[0])*(-E13 + F13 + 2*E11*d + 2*F11*d - 4*U*c*d + A11*e - E11*e - (E12*(d - e))/EE - B11*F + F11*F + 4*U*e*F - (F12*(d + F))/FF)  + 0.5*(c1[0] - c3[0] )* (A13 + E13 + B11*b - F11*b + A11*d - E11*d - 2*U*b*e - (A12*e)/AA + (E12*(d - e))/EE + 2*U*d*F - 2*U*(b*e - d*F) + (F12*b)/FF) );
//
//        fc2[1] = - invJJ*( (c2[1] - c1[1] )*( A14+E14-E11*b-F11*b+2*U*b*c+(A12*c)/(2*AA)+(E12*(b+c-2*F))/(2*EE)-A11*F+E11*F-2*U*F*F+(F12*b)/(2*FF) ) + 0.5*(c4[1] - c3[1])*(-E13 + F13 + 2*E11*d + 2*F11*d - 4*U*c*d + A11*e - E11*e - (E12*(d - e))/EE - B11*F + F11*F + 4*U*e*F - (F12*(d + F))/FF)  + 0.5*(c1[1] - c3[1] )* (A13 + E13 + B11*b - F11*b + A11*d - E11*d - 2*U*b*e - (A12*e)/AA + (E12*(d - e))/EE + 2*U*d*F - 2*U*(b*e - d*F) + (F12*b)/FF) );
//
//        fc2[2] = - invJJ*( (c2[2] - c1[2] )*( A14+E14-E11*b-F11*b+2*U*b*c+(A12*c)/(2*AA)+(E12*(b+c-2*F))/(2*EE)-A11*F+E11*F-2*U*F*F+(F12*b)/(2*FF) ) + 0.5*(c4[2] - c3[2])*(-E13 + F13 + 2*E11*d + 2*F11*d - 4*U*c*d + A11*e - E11*e - (E12*(d - e))/EE - B11*F + F11*F + 4*U*e*F - (F12*(d + F))/FF)  + 0.5*(c1[2] - c3[2] )* (A13 + E13 + B11*b - F11*b + A11*d - E11*d - 2*U*b*e - (A12*e)/AA + (E12*(d - e))/EE + 2*U*d*F - 2*U*(b*e - d*F) + (F12*b)/FF) );
//
//        fc3[0] =  - 0.5*invJJ*( (c2[0] - c1[0] )*(-A13 - F13 - B11*b + F11*b - A11*d - E11*d - 2*F11*d + 4*U*c*d - A11*e + E11*e + 2*U*b*e + (A12*e)/AA + B11*F - F11*F - 2*U*d*F - 4*U*e*F + 2*U*(b*e - d*F) - (F12*b)/FF + (F12*(d + F))/FF) + (c4[0] - c3[0] )*(-B13 - F13 + A11*a - E11*a + B11*d - F11*d - 2*U*d*e + (E12*a)/EE + 2*U*a*F - 2*U*(d*e - a*F) + (B12*F)/BB + (F12*(d + F))/FF - 2*(B14 + F14 - E11*a - F11*a + 2*U*a*c + B11*e - F11*e - 2*U*e*e + (E12*a)/(2*EE) + (B12*c)/(2*BB) + (F12*(a + c + 2*e))/(2*FF))) + (c1[0] - c3[0] ) * (-B13 - F13 + A11*a - E11*a + B11*d - F11*d - 2*U*d*e + (E12*a)/EE + 2*U*a*F - 2*U*(d*e - a*F) + (B12*F)/BB + (F12*(d + F))/FF - 2*(-2*U*((-a)*b + d*d) + (A12*a)/(2*AA) + (E12* a)/(2*EE) + (B12*b)/(2*BB) + (F12*b)/(2*FF))) );
//
//        fc3[1] =  - 0.5*invJJ*( (c2[1] - c1[1] )*(-A13 - F13 - B11*b + F11*b - A11*d - E11*d - 2*F11*d + 4*U*c*d - A11*e + E11*e + 2*U*b*e + (A12*e)/AA + B11*F - F11*F - 2*U*d*F - 4*U*e*F + 2*U*(b*e - d*F) - (F12*b)/FF + (F12*(d + F))/FF) + (c4[1] - c3[1] )*(-B13 - F13 + A11*a - E11*a + B11*d - F11*d - 2*U*d*e + (E12*a)/EE + 2*U*a*F - 2*U*(d*e - a*F) + (B12*F)/BB + (F12*(d + F))/FF - 2*(B14 + F14 - E11*a - F11*a + 2*U*a*c + B11*e - F11*e - 2*U*e*e + (E12*a)/(2*EE) + (B12*c)/(2*BB) + (F12*(a + c + 2*e))/(2*FF))) + (c1[1] - c3[1] ) * (-B13 - F13 + A11*a - E11*a + B11*d - F11*d - 2*U*d*e + (E12*a)/EE + 2*U*a*F - 2*U*(d*e - a*F) + (B12*F)/BB + (F12*(d + F))/FF - 2*(-2*U*((-a)*b + d*d) + (A12*a)/(2*AA) + (E12* a)/(2*EE) + (B12*b)/(2*BB) + (F12*b)/(2*FF))) ) ;
//
//        fc3[2] =  - 0.5*invJJ*( (c2[2] - c1[2] )*(-A13 - F13 - B11*b + F11*b - A11*d - E11*d - 2*F11*d + 4*U*c*d - A11*e + E11*e + 2*U*b*e + (A12*e)/AA + B11*F - F11*F - 2*U*d*F - 4*U*e*F + 2*U*(b*e - d*F) - (F12*b)/FF + (F12*(d + F))/FF) + (c4[2] - c3[2] )*(-B13 - F13 + A11*a - E11*a + B11*d - F11*d - 2*U*d*e + (E12*a)/EE + 2*U*a*F - 2*U*(d*e - a*F) + (B12*F)/BB + (F12*(d + F))/FF - 2*(B14 + F14 - E11*a - F11*a + 2*U*a*c + B11*e - F11*e - 2*U*e*e + (E12*a)/(2*EE) + (B12*c)/(2*BB) + (F12*(a + c + 2*e))/(2*FF))) + (c1[2] - c3[2] ) * (-B13 - F13 + A11*a - E11*a + B11*d - F11*d - 2*U*d*e + (E12*a)/EE + 2*U*a*F - 2*U*(d*e - a*F) + (B12*F)/BB + (F12*(d + F))/FF - 2*(-2*U*((-a)*b + d*d) + (A12*a)/(2*AA) + (E12* a)/(2*EE) + (B12*b)/(2*BB) + (F12*b)/(2*FF))) );
//
//
//        fc4[0] =  - invJJ*( 0.5*(c2[0] - c1[0] )*( -E13 + F13 + 2*E11*d + 2*F11*d - 4*U*c*d + A11*e - E11*e - (E12*(d - e))/EE - B11*F + F11*F +4*U*e*F - (F12*(d + F))/FF ) + (c4[0] - c3[0])*(B14 + F14 - E11*a - F11*a + 2*U*a*c + B11*e - F11*e - 2*U*e*e + (E12*a)/(2*EE) + (B12*c)/(2*BB) + (F12*(a + c + 2*e))/(2*FF))  + 0.5*(c1[0] - c3[0] )* (B13 + F13 - A11*a + E11*a - B11*d + F11*d + 2*U*d*e - (E12*a)/EE - 2*U*a*F + 2*U*(d*e - a*F) - (B12*F)/BB - (F12*(d + F))/FF) )  ;
//
//        fc4[1] =  - invJJ*( 0.5*(c2[1] - c1[1] )*( -E13 + F13 + 2*E11*d + 2*F11*d - 4*U*c*d + A11*e - E11*e - (E12*(d - e))/EE - B11*F + F11*F +4*U*e*F - (F12*(d + F))/FF ) + (c4[1] - c3[1])*(B14 + F14 - E11*a - F11*a + 2*U*a*c + B11*e - F11*e - 2*U*e*e + (E12*a)/(2*EE) + (B12*c)/(2*BB) + (F12*(a + c + 2*e))/(2*FF))  + 0.5*(c1[1] - c3[1] )* (B13 + F13 - A11*a + E11*a - B11*d + F11*d + 2*U*d*e - (E12*a)/EE - 2*U*a*F + 2*U*(d*e - a*F) - (B12*F)/BB - (F12*(d + F))/FF) ) ;
//
//        fc4[2] =  - invJJ*( 0.5*(c2[2] - c1[2] )*( -E13 + F13 + 2*E11*d + 2*F11*d - 4*U*c*d + A11*e - E11*e - (E12*(d - e))/EE - B11*F + F11*F +4*U*e*F - (F12*(d + F))/FF ) + (c4[2] - c3[2])*(B14 + F14 - E11*a - F11*a + 2*U*a*c + B11*e - F11*e - 2*U*e*e + (E12*a)/(2*EE) + (B12*c)/(2*BB) + (F12*(a + c + 2*e))/(2*FF))  + 0.5*(c1[2] - c3[2] )* (B13 + F13 - A11*a + E11*a - B11*d + F11*d + 2*U*d*e - (E12*a)/EE - 2*U*a*F + 2*U*(d*e - a*F) - (B12*F)/BB - (F12*(d + F))/FF) ) ;



//        std::cout<<a<<" "<<" "<<b<<" "<<c<<" "<<d<<" "<<e<<" "<<F<<" "<<AA<<" "<<BB<<" "<<CC<<" "<<DD<<" "<<EE<<" "
//                ""<<FF<<" "<<GG<<" "<<HH<<" "<<JJ<<" "<<ATG1<<" "<<ATG2<<" "<<ATG3<<" "<<ATG4<<" "<<U<<" "<<A1<<" "
//                ""<<A2<<" "<<E1<<" "<<E2<<" "<<B1<<" "<<B2<<" "<<F1<<" "<<F2<<" "<<A11<<" "<<A12<<" "<<A13<<" "
//                ""<<A14<<" "<<E11<<" "<<E12<<" "<<E13<<" "<<E14<<" "<<B11<<" "<<B12<<" "<<B13<<" "<<B14<<" "<<F11<<" "
//                ""<<F12<<" "<<F13<<" "<<F14<<" ";
//
//                std::cout<<fc1[0]<<" "<<fc1[1]<<" "<<fc1[2]<<" "<<fc2[0]<<" "<<fc2[1]<<" "<<fc2[2]<<" "<<fc3[0]<<" "
//                ""<<fc3[1]<<" "<<fc3[2]<<" "<<fc4[0]<<" "<<fc4[1]<<" "<<fc4[2]<<endl;
        //        " "<<ATG2<<" "
//                <<ATG3<<" "<<ATG4<<" "<<U<<endl;
#ifdef DETAILEDOUTPUT
        std::cout<<f1[0]<<" "<<f1[1]<<" "<<f1[2]<<" "<<
                 f2[0]<<" "<<f2[1]<<" "<<f2[2]<<" "<<
                 f3[0]<<" "<<f3[1]<<" "<<f3[2]<<" "<<
                 f4[0]<<" "<<f4[1]<<" "<<f4[2]<<endl;
#endif
    }
//    delete c1;
//    delete c2;
//    delete c3;
//    delete c4;
//    delete f1;
//    delete f2;
//    delete f3;
//    delete f4;
//    stop = clock();
//    elapsedtime = ((float)stop - (float)start)/CLOCKS_PER_SEC *1000;
//    std::cout<<"S CFE "<<elapsedtime<<endl;
//    CUDAvars cvars=CUDAcommon::getCUDAvars();
//    cvars.Scforce += elapsedtime;
//    CUDAcommon::cudavars=cvars;
    delete [] newc2;
}
