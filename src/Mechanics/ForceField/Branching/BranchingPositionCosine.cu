
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
#include <math.h>

#include "BranchingPositionCosine.h"
#include "BranchingPositionCosineCUDA.h"
#include "BranchingPosition.h"

#include "BranchingPoint.h"
#include "MathFunctions.h"


using namespace mathfunc;
#ifdef CUDAACCL
void BranchingPositionCosine::deallocate(){
    CUDAcommon::handleerror(cudaStreamDestroy(stream));
    CUDAcommon::handleerror(cudaFree(gU_i));
    CUDAcommon::handleerror(cudaFree(gU_sum));
}
void BranchingPositionCosine::optimalblocksnthreads( int nint){
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
                                                       BranchingPositionCosineenergy, blockToSmemFB, 0);
        blocksnthreadse.push_back((nint + blockSize - 1) / blockSize);
        blocksnthreadse.push_back(blockSize);
        blockSize = 0;

        cudaOccupancyMaxPotentialBlockSizeVariableSMem(&minGridSize, &blockSize,
                                                       BranchingPositionCosineenergyz, blockToSmemFB2, 0);
        blocksnthreadsez.push_back((nint + blockSize - 1) / blockSize);
        blocksnthreadsez.push_back(blockSize);
        blockSize = 0;

        cudaOccupancyMaxPotentialBlockSizeVariableSMem(&minGridSize, &blockSize,
                                                       BranchingPositionCosineforces, blockToSmemFB, 0);
        blocksnthreadsf.push_back((nint + blockSize - 1) / blockSize);
        blocksnthreadsf.push_back(blockSize);

        CUDAcommon::handleerror(cudaMalloc((void **) &gU_i, nint*sizeof(double)));
        CUDAcommon::handleerror(cudaMalloc((void **) &gU_sum, sizeof(double)));

        char a[] = "BranchingFF";
        char b[] = "Branching Position Cosine";
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
double* BranchingPositionCosine::energy(double *coord, double *f, int *beadSet,
                                        double *kpos, double *pos, int *params) {
    if(blocksnthreadse[1]>0) {
        BranchingPositionCosineenergy<<<blocksnthreadse[0], blocksnthreadse[1], (9 * blocksnthreadse[1]) * sizeof
                (double), stream>>> (coord, f, beadSet, kpos, pos, params, gU_i, CUDAcommon::getCUDAvars().gculpritID,
                CUDAcommon::getCUDAvars().gculpritFF,
                CUDAcommon::getCUDAvars().gculpritinteraction, gFF, ginteraction);
        auto cvars = CUDAcommon::getCUDAvars();
        cvars.streamvec.push_back(&stream);
        CUDAcommon::cudavars = cvars;
        CUDAcommon::handleerror( cudaGetLastError() ,"BranchingPositionCosineenergy", "BranchingPositionCosine.cu");
        double* gpu_Utot = CUDAcommon::getCUDAvars().gpu_energy;
        addvector<<<1,1,0,stream>>>(gU_i,params, gU_sum, gpu_Utot);
        CUDAcommon::handleerror( cudaGetLastError() ,"BranchingPositionCosineenergy", "BranchingPositionCosine.cu");
        return gU_sum;}
    else
        return NULL;
}


double* BranchingPositionCosine::energy(double *coord, double *f, int *beadSet,
                                        double *kpos, double *pos, double *z, int *params) {
    if(blocksnthreadsez[1]>0) {
        BranchingPositionCosineenergyz << < blocksnthreadsez[0], blocksnthreadsez[1], (18 * blocksnthreadsez[1]) *
                                          sizeof(double), stream>> > (coord, f, beadSet, kpos, pos, params, gU_i, z,
                CUDAcommon::getCUDAvars().gculpritID,
                CUDAcommon::getCUDAvars().gculpritFF,
                CUDAcommon::getCUDAvars().gculpritinteraction, gFF, ginteraction );
        auto cvars = CUDAcommon::getCUDAvars();
        cvars.streamvec.push_back(&stream);
        CUDAcommon::cudavars = cvars;
        CUDAcommon::handleerror(cudaGetLastError(),"BranchingPositionCosineenergyz", "BranchingPositionCosine.cu");
        double* gpu_Utot = CUDAcommon::getCUDAvars().gpu_energy;

        addvector<<<1,1,0,stream>>>(gU_i,params, gU_sum, gpu_Utot);
        CUDAcommon::handleerror(cudaGetLastError(),"BranchingPositionCosineenergyz", "BranchingPositionCosine.cu");
        return gU_sum;
    }else
        return NULL;
}

void BranchingPositionCosine::forces(double *coord, double *f, int *beadSet,
                                     double *kpos, double *pos, int *params){
    if(blocksnthreadsf[1]>0) {
        BranchingPositionCosineforces << < blocksnthreadsf[0], blocksnthreadsf[1], (9 * blocksnthreadsf[1]) *
                                                                                   sizeof(double), stream >> > (coord, f, beadSet, kpos, pos, params);
        auto cvars = CUDAcommon::getCUDAvars();
        cvars.streamvec.push_back(&stream);
        CUDAcommon::cudavars = cvars;
        CUDAcommon::handleerror(cudaGetLastError(),"BranchingPositionCosineforces", "BranchingPositionCosine.cu");
    }
}
void BranchingPositionCosine::checkforculprit() {
    CUDAcommon::printculprit("BranchingPosition","BranchingPositionCosine");
    BranchingPoint* br;
    br = (BranchingPoint::getBranchingPoints()[CUDAcommon::getCUDAvars().culpritID[0]]);
    cout<<"Printing culprit branching point information."<<endl;
    br->printSelf();
    exit(EXIT_FAILURE);
}
#endif

double BranchingPositionCosine::energy(double *coord, double *f, int *beadSet,
                                       double *kpos, double *pos){


    int n = BranchingPosition<BranchingPositionCosine>::n;
    int nint = BranchingPoint::getBranchingPoints().size();

    double *coord1, *coord2, *coord3, X, D, XD, xd, theta, posheta, dTheta, U_i;
    double *mp = new double[3];

    double U = 0;

    for(int i = 0; i < nint; i += 1) {

        coord1 = &coord[3 * beadSet[n * i]];
        coord2 = &coord[3 * beadSet[n * i + 1]];
        coord3 = &coord[3 * beadSet[n * i + 2]];

        midPointCoordinate(mp, coord1, coord2, pos[i]);
        X = sqrt(scalarProduct(mp, coord2, mp, coord2));
        D = sqrt(scalarProduct(mp, coord3, mp, coord3));

        XD = X * D;

        xd = scalarProduct(mp, coord2, mp, coord3);

        theta = safeacos(xd / XD);
        posheta = 0.5*M_PI;
        dTheta = theta-posheta;

        U_i = kpos[i] * ( 1 - cos(dTheta) );


        if(fabs(U_i) == numeric_limits<double>::infinity()
           || U_i != U_i || U_i < -1.0) {

            //set culprit and return
            BranchingInteractions::_branchingCulprit = BranchingPoint::getBranchingPoints()[i];

            return -1;
        }

        U += U_i;
    }
    delete mp;
    return U;
}

double BranchingPositionCosine::energy(double *coord, double *f, int *beadSet,
                                       double *kpos, double *pos, double d){

    int n = BranchingPosition<BranchingPositionCosine>::n;
    int nint = BranchingPoint::getBranchingPoints().size();

    double *coord1, *coord2, *coord3, *f1, *f2, *f3, X, D, XD, xd, theta, posheta, dTheta, U_i;
    double *mp = new double[3];
    double *vzero = new double[3]; vzero[0] = 0.0; vzero[1] = 0.0; vzero[2] = 0.0;

    double U = 0;

    for(int i = 0; i < nint; i += 1) {

        coord1 = &coord[3 * beadSet[n * i]];
        coord2 = &coord[3 * beadSet[n * i + 1]];
        coord3 = &coord[3 * beadSet[n * i + 2]];
        f1 = &f[3 * beadSet[n * i]];
        f2 = &f[3 * beadSet[n * i + 1]];
        f3 = &f[3 * beadSet[n * i + 2]];


        midPointCoordinateStretched(mp, coord1, f1, coord2, f2, pos[i], d);
        X = sqrt(scalarProductStretched(mp, vzero, coord2, f2, mp, vzero, coord2, f2, d));
        D = sqrt(scalarProductStretched(mp, vzero, coord3, f3, mp, vzero, coord3, f3, d));

        XD = X * D;

        xd = scalarProductStretched(mp, vzero, coord2, f2, mp, vzero, coord3, f3, d);

        theta = safeacos(xd / XD);
        posheta = 0.5*M_PI;
        dTheta = theta-posheta;

        U_i = kpos[i] * ( 1 - cos(dTheta) );
//    std::cout<<i << U_i<<endl;

        if(fabs(U_i) == numeric_limits<double>::infinity()
           || U_i != U_i || U_i < -1.0) {

            //set culprit and return
            BranchingInteractions::_branchingCulprit = BranchingPoint::getBranchingPoints()[i];

            return -1;
        }

        U += U_i;
    }
    delete mp;
    delete vzero;
    return U;
}

void BranchingPositionCosine::forces(double *coord, double *f, int *beadSet,
                                     double *kpos, double *pos){

    int n = BranchingPosition<BranchingPositionCosine>::n;
    int nint = BranchingPoint::getBranchingPoints().size();

    double *coord1, *coord2, *coord3, *f1, *f2, *f3, X, D, XD, xd, invX, invD, position, A, B, C, k, theta, posheta, dTheta, U_i;
    double *mp = new double[3];


    for(int i = 0; i < nint; i += 1) {

        coord1 = &coord[3 * beadSet[n * i]];
        coord2 = &coord[3 * beadSet[n * i + 1]];
        coord3 = &coord[3 * beadSet[n * i + 2]];
        f1 = &f[3 * beadSet[n * i]];
        f2 = &f[3 * beadSet[n * i + 1]];
        f3 = &f[3 * beadSet[n * i + 2]];

        midPointCoordinate(mp, coord1, coord2, pos[i]);
        X = sqrt(scalarProduct(mp, coord2, mp, coord2));
        D = sqrt(scalarProduct(mp, coord3, mp, coord3));

        XD = X * D;
        xd = scalarProduct(mp, coord2, mp, coord3);
//        std::cout<<xd<<" "<<scalarProduct(mp, coord2, mp, coord3)<<" "<<mp[0]<<" "<<mp[1]<<" "<<mp[2]<<" "
//                ""<<coord2[0]<<" "
//                ""<<coord2[1]<<" "<<coord2[2]<<" "
//                ""<<coord3[0]<<" "
//                ""<<coord3[1]<<" "<<coord3[2]<<" "<<endl;
        invX = 1/X;
        invD = 1/D;
        A = invX*invD;
        B = invX*invX;
        C = invD*invD;

        theta = safeacos(xd / XD);
        posheta = 0.5*M_PI;
        dTheta = theta-posheta;

        position = pos[i];

        k =  kpos[i] * A * sin(dTheta)/sin(theta);
        //bead 1
        f1[0] +=  k * (1-position)* (- (1-position)*(coord2[0] - coord1[0]) - (coord3[0] - (1-position)*coord1[0] - position*coord2[0])
                                     + xd *(B*(1-position)*(coord2[0] - coord1[0]) + C*(coord3[0] - (1-position)*coord1[0] - position*coord2[0])) );

        f1[1] +=  k * (1-position)* (- (1-position)*(coord2[1] - coord1[1]) - (coord3[1] - (1-position)*coord1[1] - position*coord2[1])
                                     + xd *(B*(1-position)*(coord2[1] - coord1[1]) + C*(coord3[1] - (1-position)*coord1[1] - position*coord2[1])) );

        f1[2] +=  k * (1-position)* (- (1-position)*(coord2[2] - coord1[2]) - (coord3[2] - (1-position)*coord1[2] - position*coord2[2])
                                     + xd *(B*(1-position)*(coord2[2] - coord1[2]) + C*(coord3[2] - (1-position)*coord1[2] - position*coord2[2])) );

        //bead 2

        f2[0] +=  k * (- position*(1-position)*(coord2[0] - coord1[0]) + (1-position)*(coord3[0]- (1-position)*coord1[0] - position*coord2[0])
                       + xd *( (1-position)*B*(1-position)*(coord2[0] - coord1[0]) - position*C*(coord3[0] - (1-position)*coord1[0] - position*coord2[0])) );

        f2[1] +=  k * (- position*(1-position)*(coord2[1] - coord1[1]) + (1-position)*(coord3[1]- (1-position)*coord1[1] - position*coord2[1])
                       + xd *( (1-position)*B*(1-position)*(coord2[1] - coord1[1]) - position*C*(coord3[1] - (1-position)*coord1[1] - position*coord2[1])) );

        f2[2] +=  k * (- position*(1-position)*(coord2[2] - coord1[2]) + (1-position)*(coord3[2]- (1-position)*coord1[2] - position*coord2[2])
                       + xd *( (1-position)*B*(1-position)*(coord2[2] - coord1[2]) - position*C*(coord3[2] - (1-position)*coord1[2] - position*coord2[2])) );

        //bead3

        f3[0] +=  k * ( (1-position)*(coord2[0] - coord1[0]) - xd * C*(coord3[0] - (1-position)*coord1[0] - position*coord2[0]) );
        f3[1] +=  k * ( (1-position)*(coord2[1] - coord1[1]) - xd * C*(coord3[1] - (1-position)*coord1[1] - position*coord2[1]) );
        f3[2] +=  k * ( (1-position)*(coord2[2] - coord1[2]) - xd * C*(coord3[2] - (1-position)*coord1[2] - position*coord2[2]) );

    }
    delete mp;
}
