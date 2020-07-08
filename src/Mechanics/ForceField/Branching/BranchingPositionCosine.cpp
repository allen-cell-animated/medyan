
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v4.0
//
//  Copyright (2015-2018)  Papoian Lab, University of Maryland
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
#include "Cylinder.h"

using namespace mathfunc;
#ifdef CUDAACCL
void BranchingPositionCosine::deallocate(){
    if(!(CUDAcommon::getCUDAvars().conservestreams))
        CUDAcommon::handleerror(cudaStreamDestroy(stream));
    CUDAcommon::handleerror(cudaFree(gU_i));
    CUDAcommon::handleerror(cudaFree(gU_sum));
    CUDAcommon::handleerror(cudaFree(gFF));
    CUDAcommon::handleerror(cudaFree(ginteraction));
}
void BranchingPositionCosine::optimalblocksnthreads( int nint){
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
        //get addition vars
        bntaddvec2.clear();
        bntaddvec2 = getaddred2bnt(nint);
        CUDAcommon::handleerror(cudaMalloc((void **) &gU_i, bntaddvec2.at(0)*sizeof(floatingpoint)));
        CUDAcommon::handleerror(cudaMemset(gU_i, 0, bntaddvec2.at(0) * sizeof(floatingpoint)));
//        CUDAcommon::handleerror(cudaMalloc((void **) &gU_i, nint*sizeof(floatingpoint)));
        CUDAcommon::handleerror(cudaMalloc((void **) &gU_sum, sizeof(floatingpoint)));

//        CUDAcommon::handleerror(cudaMalloc((void **) &gU_i, nint*sizeof(floatingpoint)));
//        CUDAcommon::handleerror(cudaMalloc((void **) &gU_sum, sizeof(floatingpoint)));

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
floatingpoint* BranchingPositionCosine::energy(floatingpoint *coord, floatingpoint *f, int *beadSet,
                                        floatingpoint *kpos, floatingpoint *pos, int *params) {
//    if(blocksnthreadse[1]>0) {
//        BranchingPositionCosineenergy<<<blocksnthreadse[0], blocksnthreadse[1], (9 * blocksnthreadse[1]) * sizeof
//                (floatingpoint), stream>>> (coord, f, beadSet, kpos, pos, params, gU_i, CUDAcommon::getCUDAvars().gculpritID,
//                CUDAcommon::getCUDAvars().gculpritFF,
//                CUDAcommon::getCUDAvars().gculpritinteraction, gFF, ginteraction);
//        auto cvars = CUDAcommon::getCUDAvars();
//        cvars.streamvec.push_back(&stream);
//        CUDAcommon::cudavars = cvars;
//        CUDAcommon::handleerror( cudaGetLastError() ,"BranchingPositionCosineenergy", "BranchingPositionCosine.cu");
//        floatingpoint* gpu_Utot = CUDAcommon::getCUDAvars().gpu_energy;
//        addvector<<<1,1,0,stream>>>(gU_i,params, gU_sum, gpu_Utot);
//        CUDAcommon::handleerror( cudaGetLastError() ,"BranchingPositionCosineenergy", "BranchingPositionCosine.cu");
//        return gU_sum;}
//    else
//        return NULL;
}


floatingpoint* BranchingPositionCosine::energy(floatingpoint *coord, floatingpoint *f, int *beadSet,
                                        floatingpoint *kpos, floatingpoint *pos, floatingpoint *z, int *params) {
    if(blocksnthreadse[1]>0) {
        BranchingPositionCosineenergy<<<blocksnthreadse[0], blocksnthreadse[1], (9 * blocksnthreadse[1]) * sizeof
                (floatingpoint), stream>>> (coord, f, beadSet, kpos, pos, params, gU_i, z, CUDAcommon::getCUDAvars()
                .gculpritID,
                CUDAcommon::getCUDAvars().gculpritFF,
                CUDAcommon::getCUDAvars().gculpritinteraction, gFF, ginteraction);
//        auto cvars = CUDAcommon::getCUDAvars();
//        cvars.streamvec.push_back(&stream);
//        CUDAcommon::cudavars = cvars;
//        CUDAcommon::handleerror( cudaGetLastError() ,"BranchingPositionCosineenergy", "BranchingPositionCosine.cu");
//        floatingpoint* gpu_Utot = CUDAcommon::getCUDAvars().gpu_energy;
//        addvector<<<1,1,0,stream>>>(gU_i,params, gU_sum, gpu_Utot);
//        CUDAcommon::handleerror( cudaGetLastError() ,"BranchingPositionCosineenergy", "BranchingPositionCosine.cu");
//        return gU_sum;
    }

    if(blocksnthreadsez[1]>0) {
        BranchingPositionCosineenergyz << < blocksnthreadsez[0], blocksnthreadsez[1], (18 * blocksnthreadsez[1]) *
                                          sizeof(floatingpoint), stream>> > (coord, f, beadSet, kpos, pos, params, gU_i, z,
                CUDAcommon::getCUDAvars().gculpritID,
                CUDAcommon::getCUDAvars().gculpritFF,
                CUDAcommon::getCUDAvars().gculpritinteraction, gFF, ginteraction );
//        auto cvars = CUDAcommon::getCUDAvars();
//        cvars.streamvec.push_back(&stream);
//        CUDAcommon::cudavars = cvars;
//        CUDAcommon::handleerror(cudaGetLastError(),"BranchingPositionCosineenergyz", "BranchingPositionCosine.cu");
//        floatingpoint* gpu_Utot = CUDAcommon::getCUDAvars().gpu_energy;
//
//        addvector<<<1,1,0,stream>>>(gU_i,params, gU_sum, gpu_Utot);
//        CUDAcommon::handleerror(cudaGetLastError(),"BranchingPositionCosineenergyz", "BranchingPositionCosine.cu");
//        return gU_sum;
    }
    if(blocksnthreadse[1]<=0 && blocksnthreadsez[1]<=0)
        return NULL;
    else{
        auto cvars = CUDAcommon::getCUDAvars();
        cvars.streamvec.push_back(&stream);
        CUDAcommon::cudavars = cvars;
        floatingpoint* gpu_Utot = CUDAcommon::getCUDAvars().gpu_energy;

//        addvector<<<1,1,0,stream>>>(gU_i,params, gU_sum, gpu_Utot);
//        cudaStreamSynchronize(stream);
//        addvectorred<<<1,200,200*sizeof(floatingpoint),stream>>>(gU_i,params, gU_sum, gpu_Utot);
//        cudaStreamSynchronize(stream);
//        std::cout<<"bntaddvec "<<bntaddvec2.at(0)<<" "<<bntaddvec2.at(1)<<" "<<bntaddvec2.at(0)<<" "
//                ""<<bntaddvec2.at(2)<<" "<<bntaddvec2.at(3)<<endl;
        resetfloatingpointvariableCUDA<<<1,1,0,stream>>>(gU_sum);
        addvectorred2<<<bntaddvec2.at(2),bntaddvec2.at(3), bntaddvec2.at(3) * sizeof(floatingpoint),stream>>>(gU_i,
                params, gU_sum, gpu_Utot);
//        CUDAcommon::handleerror(cudaDeviceSynchronize(),"FilamentBendingCosineenergyz", "FilamentBendingCosine.cu");
        CUDAcommon::handleerror(cudaGetLastError(),"FilamentBendingCosineenergyz", "FilamentBendingCosine.cu");
        return gU_sum;
    }

}

void BranchingPositionCosine::forces(floatingpoint *coord, floatingpoint *f, int *beadSet,
                                     floatingpoint *kpos, floatingpoint *pos, int *params){
    if(blocksnthreadsf[1]>0) {
        BranchingPositionCosineforces << < blocksnthreadsf[0], blocksnthreadsf[1], (9 * blocksnthreadsf[1]) *
                                                                                   sizeof(floatingpoint), stream >> > (coord, f, beadSet, kpos, pos, params);
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

floatingpoint BranchingPositionCosine::energy(floatingpoint *coord, int *beadSet,
                                       floatingpoint *kpos, floatingpoint *pos){


    int n = BranchingPosition<BranchingPositionCosine>::n;
    int nint = BranchingPoint::getBranchingPoints().size();

    floatingpoint *coord1, *coord2, *coord3, X, D, XD, xd, posheta,
    position, U_i;
    floatingpoint *mp = new floatingpoint[3];
    floatingpoint *coord2prime = new floatingpoint[3];
    floatingpoint U = 0.0;

    for(int i = 0; i < nint; i += 1) {

        coord1 = &coord[3 * beadSet[n * i]];
        coord2 = &coord[3 * beadSet[n * i + 1]];
        //If the branching point is at the plus end of the cylinder, we end up with
        // singularities in energy and force expressions. To avoid it, we will create a
        // virtual plus end and use that to define vectors.
        position = pos[i];
        if(position == 1){
            position = 0.5;//assign a dummy position and extend plus end to new coordinates.
            for(int dim = 0; dim < 3; dim++) {
                coord2prime[dim] = (1 / position) * (coord2[dim] - (1 - position) *
                                    coord1[dim]);//extended plus end coordinate;
            }
            coord2 = &coord2prime[0];
        }
        coord3 = &coord[3 * beadSet[n * i + 2]];
        midPointCoordinate(mp, coord1, coord2, position);
        X = sqrt(scalarProduct(mp, coord2, mp, coord2));
        D = sqrt(scalarProduct(mp, coord3, mp, coord3));

        XD = X * D;

        xd = scalarProduct(mp, coord2, mp, coord3);

        floatingpoint x = xd/XD;

        if(abs(abs(x) - 1.0)<0.001) {
            xd = 0.999 * XD;
            x = xd / XD;
        }

        if (x < -1.0) x = -1.0;
        else if (x > 1.0) x = 1.0;

        floatingpoint cosp =  x;
        posheta = 0.5*M_PI;
        floatingpoint sinp = sqrt(max<floatingpoint>((1-cosp*cosp),(floatingpoint)0.0));
        floatingpoint cospminusq = cosp * cos(posheta) + sinp * sin(posheta);
        U_i = kpos[i] * ( 1 - cospminusq );

        /*theta = safeacos(xd / XD);
        posheta = 0.5*M_PI;
        dTheta = theta-posheta;

        U_i = kpos[i] * ( 1 - cos(dTheta) );*/


        if(fabs(U_i) == numeric_limits<floatingpoint>::infinity()
           || U_i != U_i || U_i < -1.0) {

            //set culprit and return
            BranchingInteractions::_branchingCulprit = BranchingPoint::getBranchingPoints()[i];

            return -1;
        }

        U += U_i;
    }
    delete[] mp;
    delete[] coord2prime;
    return U;
}

floatingpoint BranchingPositionCosine::energy(floatingpoint *coord, floatingpoint *f, int *beadSet,
                                       floatingpoint *kpos, floatingpoint *pos, floatingpoint d){

    int n = BranchingPosition<BranchingPositionCosine>::n;
    int nint = BranchingPoint::getBranchingPoints().size();

    floatingpoint *coord1, *coord2, *coord3, X, D, XD, xd, posheta, U_i;
    floatingpoint *f1, *f2, *f3;
    floatingpoint *mp = new floatingpoint[3];
    floatingpoint *vzero = new floatingpoint[3]; vzero[0] = 0.0; vzero[1] = 0.0; vzero[2] = 0.0;

    floatingpoint U = 0.0;

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

        floatingpoint x = xd/XD;

        if(abs(abs(x) - 1.0)<0.001) {
            xd = 0.999 * XD;
            x = xd / XD;
        }

        if (x < -1.0) x = -1.0;
        else if (x > 1.0) x = 1.0;

        floatingpoint cosp =  x;
        posheta = 0.5*M_PI;
        floatingpoint sinp = sqrt(max<floatingpoint>((1-cosp*cosp),(floatingpoint)0.0));
        floatingpoint cospminusq = cosp * cos(posheta) + sinp * sin(posheta);
        U_i = kpos[i] * ( 1 - cospminusq );

        /*theta = safeacos(xd / XD);
        posheta = 0.5*M_PI;
        dTheta = theta-posheta;

        U_i = kpos[i] * ( 1 - cos(dTheta) );*/

        if(fabs(U_i) == numeric_limits<floatingpoint>::infinity()
           || U_i != U_i || U_i < -1.0) {

            //set culprit and return
            BranchingInteractions::_branchingCulprit = BranchingPoint::getBranchingPoints()[i];

            return -1;
        }

        U += U_i;
    }
    delete[] mp;
    delete[] vzero;
    return U;
}

void BranchingPositionCosine::forces(floatingpoint *coord, floatingpoint *f, int *beadSet,
                                     floatingpoint *kpos, floatingpoint *pos,
                                     floatingpoint *stretchforce){

    int n = BranchingPosition<BranchingPositionCosine>::n;
    int nint = BranchingPoint::getBranchingPoints().size();

    floatingpoint *coord1, *coord2, *coord3, X, D, XD, xd, invX, invD, position, A, B, C, k, posheta;
	floatingpoint  *f1, *f2, *f3;
    floatingpoint *mp = new floatingpoint[3];
    floatingpoint *coord2prime = new floatingpoint[3];
    floatingpoint f1tempx, f1tempy, f1tempz, f2tempx, f2tempy, f2tempz, f3tempx, f3tempy,
    f3tempz;

    for(int i = 0; i < nint; i += 1) {

        coord1 = &coord[3 * beadSet[n * i]];
        coord2 = &coord[3 * beadSet[n * i + 1]];
        //If the branching point is at the plus end of the cylinder, we end up with
        // singularities in energy and force expressions. To avoid it, we will create a
        // virtual plus end and use that to define vectors.
        position = pos[i];
        if(areEqual(position, (floatingpoint)1.0)){
            position = 0.5;//assign a dummy position and extend plus end to new coordinates.
            for(int dim = 0; dim < 3; dim++) {
                coord2prime[dim] = (1 / position) * (coord2[dim] - (1 - position) *
                                                                   coord1[dim]);//extended plus end coordinate;
            }
            coord2 = &coord2prime[0];
        }
        coord3 = &coord[3 * beadSet[n * i + 2]];
        f1 = &f[3 * beadSet[n * i]];
        f2 = &f[3 * beadSet[n * i + 1]];
        f3 = &f[3 * beadSet[n * i + 2]];

        midPointCoordinate(mp, coord1, coord2, position);
        X = sqrt(scalarProduct(mp, coord2, mp, coord2));
        D = sqrt(scalarProduct(mp, coord3, mp, coord3));

        XD = X * D;
        xd = scalarProduct(mp, coord2, mp, coord3);
        invX = 1/X;
        invD = 1/D;
        A = invX*invD;
        B = invX*invX;
        C = invD*invD;

        floatingpoint x = xd/XD;

        if(abs(abs(x) - 1.0)<0.001) {
            xd = 0.999 * XD;
            x = xd / XD;
        }

        if (x < -1.0) x = -1.0;
        else if (x > 1.0) x = 1.0;

        floatingpoint cosp =  x;
        posheta = 0.5*M_PI;
        floatingpoint sinp = sqrt(max<floatingpoint>((1-cosp*cosp),(floatingpoint)0.0));
        floatingpoint sinpminusq = sinp * cos(posheta) - cosp * sin(posheta);

        k = kpos[i] * A * sinpminusq/sinp;

	    /*if(abs(xd/XD - 1.0)<0.001){
		    xd = 0.999*XD;
	    }

        theta = safeacos(xd / XD);
        posheta = 0.5*M_PI;
        dTheta = theta-posheta;

        position = pos[i];

        k =  kpos[i] * A * sin(dTheta)/sin(theta);*/

	    //If the branching point is NOT bound to plusend, the local force variables
	    // ftemp1[3], ftemp2[3], ftemp3[3] represent forces on parent_minus, parent_plus
	    // and offspring_minus ends respectively.
	    // If the branching point IS bound to plusend, the force variables represent
	    // forces on parent_minus, parent_extendedplus and offspring_minus ends
	    // respectively. Under this condition, a transformation is necessary to realize
	    // the actual forces on the beads of interest.

        //bead 1
        f1tempx =  k * (1-position)* (- (1-position)*(coord2[0] - coord1[0]) - (coord3[0]
            - (1-position)*coord1[0] - position*coord2[0])
                                     + xd *(B*(1-position)*(coord2[0] - coord1[0]) + C*(coord3[0] - (1-position)*coord1[0] - position*coord2[0])) );

        f1tempy =  k * (1-position)* (- (1-position)*(coord2[1] - coord1[1]) - (coord3[1]
            - (1-position)*coord1[1] - position*coord2[1])
                                     + xd *(B*(1-position)*(coord2[1] - coord1[1]) + C*(coord3[1] - (1-position)*coord1[1] - position*coord2[1])) );

        f1tempz =  k * (1-position)* (- (1-position)*(coord2[2] - coord1[2]) - (coord3[2]
            - (1-position)*coord1[2] - position*coord2[2])
                                     + xd *(B*(1-position)*(coord2[2] - coord1[2]) + C*(coord3[2] - (1-position)*coord1[2] - position*coord2[2])) );

        //bead 2
        f2tempx =  k * (- position*(1-position)*(coord2[0] - coord1[0]) + (1-position)*
            (coord3[0]- (1-position)*coord1[0] - position*coord2[0])
                       + xd *( (1-position)*B*(1-position)*(coord2[0] - coord1[0]) - position*C*(coord3[0] - (1-position)*coord1[0] - position*coord2[0])) );

        f2tempy =  k * (- position*(1-position)*(coord2[1] - coord1[1]) + (1-position)*
            (coord3[1]- (1-position)*coord1[1] - position*coord2[1])
                       + xd *( (1-position)*B*(1-position)*(coord2[1] - coord1[1]) - position*C*(coord3[1] - (1-position)*coord1[1] - position*coord2[1])) );

        f2tempz =  k * (- position*(1-position)*(coord2[2] - coord1[2]) + (1-position)*
            (coord3[2]- (1-position)*coord1[2] - position*coord2[2])
                       + xd *( (1-position)*B*(1-position)*(coord2[2] - coord1[2]) - position*C*(coord3[2] - (1-position)*coord1[2] - position*coord2[2])) );

        //bead3
        f3tempx =  k * ( (1-position)*(coord2[0] - coord1[0]) - xd * C*
                (coord3[0] - (1-position)*coord1[0] - position*coord2[0]) );
        f3tempy =  k * ( (1-position)*(coord2[1] - coord1[1]) - xd * C*
                (coord3[1] - (1-position)*coord1[1] - position*coord2[1]) );
        f3tempz =  k * ( (1-position)*(coord2[2] - coord1[2]) - xd * C*
                (coord3[2] - (1-position)*coord1[2] - position*coord2[2]) );

        f3[0] += f3tempx;
        f3[1] += f3tempy;
        f3[2] += f3tempz;

        stretchforce[3*i] = f3tempx;
        stretchforce[3*i + 1] = f3tempy;
        stretchforce[3*i + 2] = f3tempz;

        //If you had calculated forces on the extended plus end, additional
        // transformations are needed.
        if(areEqual(pos[i],(floatingpoint)1.0)){
            floatingpoint factor = (position-1)/position;
            f1[0] += f1tempx + f2tempx*factor;
            f1[1] += f1tempy + f2tempy*factor;
            f1[2] += f1tempz + f2tempz*factor;

            f2[0] += f2tempx*(1/position);
            f2[1] += f2tempy*(1/position);
            f2[2] += f2tempz*(1/position);
        }
        else{
            f1[0] += f1tempx;
            f1[1] += f1tempy;
            f1[2] += f1tempz;

            f2[0] += f2tempx;
            f2[1] += f2tempy;
            f2[2] += f2tempz;
        }


	    #ifdef CHECKFORCES_INF_NAN
	    if(checkNaN_INF<floatingpoint>(f1, 0, 2)||checkNaN_INF<floatingpoint>(f2,0,2)||checkNaN_INF<floatingpoint>(f3,0,2)){
		    cout<<"Branching Position Force becomes infinite. Printing data "<<endl;

            auto b = BranchingPoint::getBranchingPoints()[i];
            auto cyl1 = b->getFirstCylinder();
            auto cyl2 = b->getSecondCylinder();
            cout<<"Cylinder IDs "<<cyl1->getId()<<" "<<cyl2->getId()<<" with cIndex "
                <<cyl1->getStableIndex()<<" "<<cyl2->getStableIndex()<<" and bIndex "
                <<cyl1->getFirstBead()->getStableIndex()<<" "
                <<cyl1->getSecondBead()->getStableIndex()<<" "
                <<cyl2->getFirstBead()->getStableIndex()<<" "
                <<cyl2->getSecondBead()->getStableIndex()<<endl;
            cyl1->adjustedrelativeposition(pos[i], true);
            cout<<"Printing intermediate variables"<<endl;
            cout<<"XD="<<XD<<", xd="<<xd<<", invX="<<invX<<", invD="<<invD<<", A="<<A
                <<", B="<<B<<", C="<<C<<", cosp="<<cosp<<", sinp="<<sinp
                <<", sinpminusq="<<sinpminusq<<", position="<<position<<endl;

		    cout<<"Printing coords"<<endl;
		    cout<<coord1[0]<<" "<<coord1[1]<<" "<<coord1[2]<<endl;
		    cout<<coord2[0]<<" "<<coord2[1]<<" "<<coord2[2]<<endl;
		    cout<<coord3[0]<<" "<<coord3[1]<<" "<<coord3[2]<<endl;

		    cout<<"Printing force"<<endl;
		    cout<<f1[0]<<" "<<f1[1]<<" "<<f1[2]<<endl;
		    cout<<f2[0]<<" "<<f2[1]<<" "<<f2[2]<<endl;
		    cout<<f3[0]<<" "<<f3[1]<<" "<<f3[2]<<endl;

		    cout<<"Printing binary Coords"<<endl;
		    printvariablebinary(coord1,0,2);
		    printvariablebinary(coord2,0,2);
		    printvariablebinary(coord3,0,2);

		    cout<<"Printing binary Force"<<endl;
		    printvariablebinary(f1,0,2);
		    printvariablebinary(f2,0,2);
		    printvariablebinary(f3,0,2);

		    exit(EXIT_FAILURE);
	    }
	    #endif
    }
    delete[] mp;
    delete[] coord2prime;
}
