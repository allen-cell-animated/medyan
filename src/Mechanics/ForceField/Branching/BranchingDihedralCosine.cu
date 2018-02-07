
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

#include "BranchingDihedralCosine.h"
#include "BranchingDihedralCosineCUDA.h"
#include "BranchingDihedral.h"

#include "BranchingPoint.h"
#include "Bead.h"

#include "MathFunctions.h"
#include <cuda.h>
#include <cuda_runtime.h>
#include "nvToolsExt.h"

using namespace mathfunc;
#ifdef CUDAACCL
void BranchingDihedralCosine::deallocate(){
    CUDAcommon::handleerror(cudaStreamDestroy(stream));
    CUDAcommon::handleerror(cudaFree(gU_i));
    CUDAcommon::handleerror(cudaFree(gU_sum));
}
void BranchingDihedralCosine::optimalblocksnthreads( int nint){
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
                                                       BranchingDihedralCosineenergy, blockToSmem, 0);
//    std::cout<<(nint +blockSize -1) / blockSize<<" "<<blockSize<<endl;
//
//    cudaOccupancyMaxPotentialBlockSize( &minGridSize, &blockSize,
//                                        CUDAExclVolRepulsionenergy, 0, 0);
        blocksnthreadse.push_back((nint + blockSize - 1) / blockSize);
        blocksnthreadse.push_back(blockSize);
//    std::cout<<(nint +blockSize -1) / blockSize<<" "<<blockSize<<endl;
        blockSize = 0;

        cudaOccupancyMaxPotentialBlockSizeVariableSMem(&minGridSize, &blockSize,
                                                       BranchingDihedralCosineenergyz, blockToSmemez, 0);
        blocksnthreadsez.push_back((nint + blockSize - 1) / blockSize);
        blocksnthreadsez.push_back(blockSize);
        blockSize = 0;

        cudaOccupancyMaxPotentialBlockSizeVariableSMem(&minGridSize, &blockSize,
                                                       BranchingDihedralCosineforces, blockToSmem, 0);
        blocksnthreadsf.push_back((nint + blockSize - 1) / blockSize);
        blocksnthreadsf.push_back(blockSize);

        CUDAcommon::handleerror(cudaMalloc((void **) &gU_i, nint*sizeof(double)));
        CUDAcommon::handleerror(cudaMalloc((void **) &gU_sum, sizeof(double)));
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
double* BranchingDihedralCosine::energy(double *coord, double *f, int *beadSet,
                                         double *kdih, double *pos, int *params) {
    if(blocksnthreadse[1]>0) {

        BranchingDihedralCosineenergy<<<blocksnthreadse[0], blocksnthreadse[1], (12 * blocksnthreadse[1]) * sizeof
                (double), stream>>>
                          (coord, f, beadSet, kdih, pos, params, gU_i);
        auto cvars = CUDAcommon::getCUDAvars();
        cvars.streamvec.push_back(&stream);
        CUDAcommon::cudavars = cvars;
        CUDAcommon::handleerror( cudaGetLastError(),"BranchingDihedralCosineenergy", "BranchingDihedralCosine.cu");
        double* gpu_Utot = CUDAcommon::getCUDAvars().gpu_energy;
        addvector<<<1,1,0,stream>>>(gU_i, params, gU_sum, gpu_Utot);
        CUDAcommon::handleerror( cudaGetLastError() ,"BranchingDihedralCosineenergy", "BranchingDihedralCosine.cu");
        return gU_sum;}
    else
        return NULL;
}


double* BranchingDihedralCosine::energy(double *coord, double *f, int *beadSet,
                                         double *kdih, double *pos, double *z,
                                         int *params) {

    if(blocksnthreadsez[1]>0) {
        BranchingDihedralCosineenergyz << < blocksnthreadsez[0], blocksnthreadsez[1], (24 * blocksnthreadsez[1]) *
                                            sizeof(double), stream>> > (coord, f, beadSet, kdih, pos,
                                            params, gU_i, z );
        auto cvars = CUDAcommon::getCUDAvars();
        cvars.streamvec.push_back(&stream);
        CUDAcommon::cudavars = cvars;
        CUDAcommon::handleerror(cudaGetLastError(),"BranchingDihedralCosineenergyz", "BranchingDihedralCosine.cu");
        double* gpu_Utot = CUDAcommon::getCUDAvars().gpu_energy;
        addvector<<<1,1,0,stream>>>(gU_i, params, gU_sum, gpu_Utot);
        CUDAcommon::handleerror(cudaGetLastError(),"BranchingDihedralCosineenergyz", "BranchingDihedralCosine.cu");

        return gU_sum;
    }else
        return NULL;
}
void BranchingDihedralCosine::forces(double *coord, double *f, int *beadSet,
                                      double *kdih,  double *pos, int *params) {
    if (blocksnthreadsf[1] > 0) {
        BranchingDihedralCosineforces << < blocksnthreadsf[0], blocksnthreadsf[1], (12 * blocksnthreadsf[1]) *
                                            sizeof(double), stream >> > (coord, f, beadSet, kdih, pos, params);
        auto cvars = CUDAcommon::getCUDAvars();
        cvars.streamvec.push_back(&stream);
        CUDAcommon::cudavars = cvars;
        CUDAcommon::handleerror(cudaGetLastError(),"BranchingDihedralCosineforces", "BranchingDihedralCosine.cu");
    }
}
#endif
double BranchingDihedralCosine::energy(double *coord, double *f, int *beadSet,
                                       double *kdih, double *pos){


    int n = BranchingDihedral<BranchingDihedralCosine>::n;
    int nint = BranchingPoint::getBranchingPoints().size();


    double *coord1, *coord2, *coord3, *coord4, n1n2, U_i;
    double *mp = new double[3];
    double *n1 = new double[3];
    double *n2 = new double[3];

    double U = 0;

    for(int i = 0; i < nint; i += 1) {

        coord1 = &coord[3 * beadSet[n * i]];
        coord2 = &coord[3 * beadSet[n * i + 1]];
        coord3 = &coord[3 * beadSet[n * i + 2]];
        coord4 = &coord[3 * beadSet[n * i + 3]];

        midPointCoordinate(mp, coord1, coord2, pos[i]);

        vectorProduct(n1, mp, coord2, mp, coord3);
        vectorProduct(n2, coord3, coord4, mp, coord3);

        normalizeVector(n1);
        normalizeVector(n2);
        n1n2 = dotProduct(n1, n2);

        U_i = kdih[i] * ( 1 - n1n2 );

        if(fabs(U_i) == numeric_limits<double>::infinity()
           || U_i != U_i || U_i < -1.0) {

            //set culprit and return
            BranchingInteractions::_branchingCulprit = BranchingPoint::getBranchingPoints()[i];

            return -1;
        }

        U += U_i;
    }
    delete mp;
    delete n1;
    delete n2;

    return U;
}


double BranchingDihedralCosine::energy(double *coord, double *f, int *beadSet,
                                       double *kdih, double *pos, double d){

    int n = BranchingDihedral<BranchingDihedralCosine>::n;
    int nint = BranchingPoint::getBranchingPoints().size();


    double *coord1, *coord2, *coord3, *coord4, *f1, *f2, *f3, *f4, n1n2, U_i;
    double *mp = new double[3];
    double *n1 = new double[3];
    double *n2 = new double[3];
    double *zero = new double[3]; zero[0] = 0; zero[1] = 0; zero[2] = 0;

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

        midPointCoordinateStretched(mp, coord1, f1, coord2, f2, pos[i], d);

        vectorProductStretched(n1, mp, zero, coord2, f2, mp, zero, coord3, f3, d);
        vectorProductStretched(n2, coord3, f3, coord4, f4, mp, zero, coord3, f3, d);

        normalizeVector(n1);
        normalizeVector(n2);
        n1n2 = dotProduct(n1, n2);

        U_i = kdih[i] * ( 1 - n1n2 );

        if(fabs(U_i) == numeric_limits<double>::infinity()
           || U_i != U_i || U_i < -1.0) {

            //set culprit and return
            BranchingInteractions::_branchingCulprit = BranchingPoint::getBranchingPoints()[i];

            return -1;
        }

        U += U_i;
    }
    delete mp;
    delete n1;
    delete n2;
    delete zero;
    return U;
}

void BranchingDihedralCosine::forces(double *coord, double *f, int *beadSet,
                                     double *kdih, double *pos){

    int n = BranchingDihedral<BranchingDihedralCosine>::n;
    int nint = BranchingPoint::getBranchingPoints().size();


    double *coord1, *coord2, *coord3, *coord4, *f1, *f2, *f3, *f4, N1, N2, n1n2, f0, NN1, NN2, X, D, Y, position;
    double n2x, n1y, xd, yd, xx, xy, yy, XD, X1, X2, Y1, Y2, D1, D2, YD;
    double *mp = new double[3];
    double *n1 = new double[3];
    double *n2 = new double[3];
    double *zero = new double[3]; zero[0] = 0; zero[1] = 0; zero[2] = 0;

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

        midPointCoordinate(mp, coord1, coord2, pos[i]);

        vectorProduct(n1, mp, coord2, mp, coord3);
        vectorProduct(n2, coord3, coord4, mp, coord3);

        N1 = sqrt(dotProduct(n1, n1));
        N2 = sqrt(dotProduct(n2, n2));
        n1n2 = dotProduct(n1, n2);

        f0 = kdih[i]/N1/N2;

        NN1 = n1n2/N1/N1;
        NN2 = n1n2/N2/N2;

        X = sqrt(scalarProduct(mp, coord2, mp, coord2));
        D = sqrt(scalarProduct(mp, coord3, mp, coord3));
        Y = sqrt(scalarProduct(coord3, coord4, coord3, coord4));

        n2x = scalarProduct(zero, n2, mp, coord2);
        n1y = scalarProduct(zero, n1, coord3, coord4);
        xd = scalarProduct(mp, coord2, mp, coord3);
        yd = scalarProduct(coord3, coord4, mp, coord3);

        xx = scalarProduct(mp, coord2, mp, coord2);
        xy = scalarProduct(mp, coord2, coord3, coord4);
        yy = scalarProduct(coord3, coord4, coord3, coord4);

        XD = n2x/D/X/X/X;
        X1 = -NN2*xd/D/X + yd/D/Y + yd/D/D/X/Y;
        X2 = xd*yd/D/D/X/X/X/Y;
        Y1 = -xd/D/X - xd/D/D/X/Y + NN1*yd/D/Y;
        Y2 = xd*yd/D/D/X/Y/Y/Y;
        D1 = NN2*xx/D/X - xy/D/X-xy/D/Y - 2*xy/D/D/X/Y + NN1*yy/D/Y;
        D2 = xd*xy/D/D/X/X/X/Y;
        YD = n1y/D/Y/Y/Y;

        position = pos[i];

        //force on b1:
        f1[0] += f0*(- (1 - position)*XD*(1-position)*( (coord2[1] - coord1[1])*(coord3[2] - (1-position)*coord1[2] - position*coord2[2]) - (coord2[2] - coord1[2])*(coord3[1] - (1-position)*coord1[1] - position*coord2[1]) ) + (1 - position)*(X1 - X2)*(1-position)*(coord2[0] - coord1[0]) - (1 - position)*Y1*(coord4[0] - coord3[0]) + (1 - position)*(D1 + D2)*(coord3[0] - (1-position)*coord1[0] - position*coord2[0]));


        f1[1] += f0*(- (1 - position)*XD*(1-position)*( (coord2[2] - coord1[2])*(coord3[0] - (1-position)*coord1[0] - position*coord2[0]) - (coord2[0] - coord1[0])*(coord3[2] - (1-position)*coord1[2] - position*coord2[2]) ) + (1 - position)*(X1 - X2)*(1-position)*(coord2[1] - coord1[1]) - (1 - position)*Y1*(coord4[1] - coord3[1]) + (1 - position)*(D1 + D2)*(coord3[1] - (1-position)*coord1[1] - position*coord2[1]));

        f1[2] += f0*(- (1 - position)*XD*(1-position)*( (coord2[0] - coord1[0])*(coord3[1] - (1-position)*coord1[1] - position*coord2[1]) - (coord2[1] - coord1[1])*(coord3[0] - (1-position)*coord1[0] - position*coord2[0]) ) + (1 - position)*(X1 - X2)*(1-position)*(coord2[2] - coord1[2]) - (1 - position)*Y1*(coord4[2] - coord3[2]) + (1 - position)*(D1 + D2)*(coord3[2] - (1-position)*coord1[2] - position*coord2[2]));


        //force on b2:
        f2[0] += f0*( (1 - position)*XD*(1-position)*( (coord2[1] - coord1[1])*(coord3[2] - (1-position)*coord1[2] - position*coord2[2]) - (coord2[2] - coord1[2])*(coord3[1] - (1-position)*coord1[1] - position*coord2[1]) ) + (X2 + position*(X1 - X2))*(1-position)*(coord2[0] - coord1[0]) - position*Y1*(coord4[0] - coord3[0]) + (position*(D1 + D2) - D2)*(coord3[0] - (1-position)*coord1[0] - position*coord2[0]) );

        f2[1] += f0*( (1 - position)*XD*(1-position)*( (coord2[2] - coord1[2])*(coord3[0] - (1-position)*coord1[0] - position*coord2[0]) - (coord2[0] - coord1[0])*(coord3[2] - (1-position)*coord1[2] - position*coord2[2]) ) + (X2 + position*(X1 - X2))*(1-position)*(coord2[1] - coord1[1]) - position*Y1*(coord4[1] - coord3[1]) + (position*(D1 + D2) - D2)*(coord3[1] - (1-position)*coord1[1] - position*coord2[1]) );

        f2[2] += f0*( (1 - position)*XD*(1-position)*( (coord2[0] - coord1[0])*(coord3[1] - (1-position)*coord1[1] - position*coord2[1]) - (coord2[1] - coord1[1])*(coord3[0] - (1-position)*coord1[0] - position*coord2[0]) ) + (X2 + position*(X1 - X2))*(1-position)*(coord2[2] - coord1[2]) - position*Y1*(coord4[2] - coord3[2]) + (position*(D1 + D2) - D2)*(coord3[2] - (1-position)*coord1[2] - position*coord2[2]) );

        //force on b3:
        f3[0] += f0*(-YD*( (coord4[1] - coord3[1])*(coord3[2] - (1-position)*coord1[2] - position*coord2[2]) - (coord4[2] - coord3[2])*(coord3[1] - (1-position)*coord1[1] - position*coord2[1]) ) - X1*(1-position)*(coord2[0] - coord1[0]) + (Y1 - Y2)*(coord4[0] - coord3[0]) + (D2 - D1)*(coord3[0] - (1-position)*coord1[0] - position*coord2[0]));

        f3[1] += f0*(-YD*( (coord4[2] - coord3[2])*(coord3[0] - (1-position)*coord1[0] - position*coord2[0]) - (coord4[0] - coord3[0])*(coord3[2] - (1-position)*coord1[2] - position*coord2[2]) ) - X1*(1-position)*(coord2[1] - coord1[1]) + (Y1 - Y2)*(coord4[1] - coord3[1]) + (D2 - D1)*(coord3[1] - (1-position)*coord1[1] - position*coord2[1]));

        f3[2] += f0*(-YD*( (coord4[0] - coord3[0])*(coord3[1] - (1-position)*coord1[1] - position*coord2[1]) - (coord4[1] - coord3[1])*(coord3[0] - (1-position)*coord1[0] - position*coord2[0]) ) - X1*(1-position)*(coord2[2] - coord1[2]) + (Y1 - Y2)*(coord4[2] - coord3[2]) + (D2 - D1)*(coord3[2] - (1-position)*coord1[2] - position*coord2[2]));


        //force on b4:
        f4[0] +=f0*( YD*( (coord4[1] - coord3[1])*(coord3[2] - (1-position)*coord1[2] - position*coord2[2]) - (coord4[2] - coord3[2])*(coord3[1] - (1-position)*coord1[1] - position*coord2[1]) ) + Y2*(coord4[0] - coord3[0]) - D2*(coord3[0] - (1-position)*coord1[0] - position*coord2[0]) );

        f4[1] +=f0*( YD*( (coord4[2] - coord3[2])*(coord3[0] - (1-position)*coord1[0] - position*coord2[0]) - (coord4[0] - coord3[0])*(coord3[2] - (1-position)*coord1[2] - position*coord2[2]) ) + Y2*(coord4[1] - coord3[1]) - D2*(coord3[1] - (1-position)*coord1[1] - position*coord2[1]) );

        f4[2] +=f0*( YD*( (coord4[0] - coord3[0])*(coord3[1] - (1-position)*coord1[1] - position*coord2[1]) - (coord4[1] - coord3[1])*(coord3[0] - (1-position)*coord1[0] - position*coord2[0]) ) + Y2*(coord4[2] - coord3[2]) - D2*(coord3[2] - (1-position)*coord1[2] - position*coord2[2]) );
    }
    delete mp;
    delete n1;
    delete n2;
    delete zero;
}
