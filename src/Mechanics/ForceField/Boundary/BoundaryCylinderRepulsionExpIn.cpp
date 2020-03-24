
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

#include "BoundaryCylinderRepulsionExpIn.h"
#include "BoundaryCylinderRepulsionIn.h"

#include "BoundaryElement.h"
#include "Bead.h"
//TODO added for temporary CUDA force.
#include "CGMethod.h"
#include "BoundaryCylinderRepulsionExpCUDA.h"

#include "cross_check.h"
#include "Cylinder.h"
#include "MathFunctions.h"
#ifdef CUDAACCL
#include "nvToolsExt.h"
#include "CUDAcommon.h"
#endif
using namespace mathfunc;
#ifdef CUDAACCL
void BoundaryCylinderRepulsionExpIn::deallocate(){
    //    if(!(CUDAcommon::getCUDAvars().conservestreams))
    //        CUDAcommon::handleerror(cudaStreamDestroy(stream));
    CUDAcommon::handleerror(cudaFree(gU_i));
    CUDAcommon::handleerror(cudaFree(gU_sum));
    CUDAcommon::handleerror(cudaFree(gFF));
    CUDAcommon::handleerror(cudaFree(ginteraction));
    //Memory alloted
    //@{
    //    size_t allocmem = 0;
    //    allocmem += (bntaddvec2.at(0) + 1) * sizeof(floatingpoint) + 200 * sizeof(char);
    //    auto c = CUDAcommon::getCUDAvars();
    //    c.memincuda -= allocmem;
    //    CUDAcommon::cudavars = c;
    //    std::cout<<"Total allocated memory "<<c.memincuda/1024<<endl;
    //    std::cout<<"Memory allocated 0 . Memory freed "<<allocmem/1024<<endl;
    //@}
}
void BoundaryCylinderRepulsionExpIn::optimalblocksnthreads( int nint, cudaStream_t
                                                         stream_pass){
    //    //CUDA stream create
    //    if(stream == NULL || !(CUDAcommon::getCUDAvars().conservestreams))
    //        CUDAcommon::handleerror(cudaStreamCreate(&stream));
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
                                                       BoundaryCylinderRepulsionExpenergy, blockToSmemF, 0);
        blocksnthreadse.push_back((nint + blockSize - 1) / blockSize);
        blocksnthreadse.push_back(blockSize);
        blockSize = 0;
        
        cudaOccupancyMaxPotentialBlockSizeVariableSMem(&minGridSize, &blockSize,
                                                       //                                                       BoundaryCylinderRepulsionExpenergyz, blockToSmemFB3, 0);
                                                       BoundaryCylinderRepulsionExpenergyz, blockToSmemZero, 0);
        blocksnthreadsez.push_back((nint + blockSize - 1) / blockSize);
        
        blocksnthreadsez.push_back(blockSize);
        blockSize = 0;
        
        cudaOccupancyMaxPotentialBlockSizeVariableSMem(&minGridSize, &blockSize,
                                                       //                                                       BoundaryCylinderRepulsionExpforces, blockToSmemF, 0);
                                                       BoundaryCylinderRepulsionExpforces, blockToSmemZero, 0);
        blocksnthreadsf.push_back((nint + blockSize - 1) / blockSize);
        blocksnthreadsf.push_back(blockSize);
        
        //get addition vars
        bntaddvec2.clear();
        bntaddvec2 = getaddred2bnt(nint);
        CUDAcommon::handleerror(cudaMalloc((void **) &gU_i, bntaddvec2.at(0)*sizeof(floatingpoint)));
        CUDAcommon::handleerror(cudaMemsetAsync(gU_i, 0, bntaddvec2.at(0) * sizeof(floatingpoint), stream));
        CUDAcommon::handleerror(cudaMalloc((void **) &gU_sum, sizeof(floatingpoint)));
        char a[] = "Boundary Cylinder Repulsion";
        char b[] = "Boundary Cylinder Repulsion Exp";
        CUDAcommon::handleerror(cudaMalloc((void **) &gFF, 100 * sizeof(char)));
        CUDAcommon::handleerror(cudaMalloc((void **) &ginteraction, 100 * sizeof(char)));
        CUDAcommon::handleerror(cudaMemcpyAsync(gFF, a, 100 * sizeof(char),
                                                cudaMemcpyHostToDevice, stream));
        CUDAcommon::handleerror(cudaMemcpyAsync(ginteraction, b, 100 * sizeof(char),
                                                cudaMemcpyHostToDevice, stream));
        //Memory alloted
        //@{
        //        size_t allocmem = 0;
        //        allocmem += (bntaddvec2.at(0) + 1) * sizeof(floatingpoint) + 200 * sizeof(char);
        //        auto c = CUDAcommon::getCUDAvars();
        //        c.memincuda += allocmem;
        //        CUDAcommon::cudavars = c;
        //        std::cout<<"Total allocated memory "<<c.memincuda/1024<<endl;
        //        std::cout<<"Memory allocated "<< allocmem/1024<<"Memory freed 0"<<endl;
        //@}
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
floatingpoint* BoundaryCylinderRepulsionExpIn::energy(floatingpoint *coord, floatingpoint *f, int *beadSet, floatingpoint *krep, floatingpoint *slen,
                                                    int* nintvec, floatingpoint* beListplane, int *params) {
    //    if(blocksnthreadse[1]>0) {
    //
    //        BoundaryCylinderRepulsionExpenergy<<<blocksnthreadse[0], blocksnthreadse[1], (3 * blocksnthreadse[1]) * sizeof
    //                (floatingpoint), stream>>>
    //                          (coord, f, beadSet, krep, slen, nintvec, beListplane, params, gU_i, CUDAcommon::getCUDAvars()
    //                          .gculpritID, CUDAcommon::getCUDAvars().gculpritFF,
    //                          CUDAcommon::getCUDAvars().gculpritinteraction, gFF, ginteraction);
    //        auto cvars = CUDAcommon::getCUDAvars();
    //        cvars.streamvec.push_back(&stream);
    //        CUDAcommon::cudavars = cvars;
    //        floatingpoint* gpu_Utot = CUDAcommon::getCUDAvars().gpu_energy;
    //        addvector<<<1,1,0,stream>>>(gU_i,params, gU_sum, gpu_Utot);
    //        CUDAcommon::handleerror( cudaGetLastError() ,"BoundaryCylinderRepulsionExpenergy", "BoundaryCylinderRepulsionExp.cu");
    //        return gU_sum;}
    //    else
    //        return NULL;
}

floatingpoint* BoundaryCylinderRepulsionExpIn::energy(floatingpoint *coord, floatingpoint *f, int *beadSet, floatingpoint *krep, floatingpoint *slen,
                                                    int* nintvec, floatingpoint* beListplane, floatingpoint *z,int *params) {
    
    //    if(blocksnthreadse[1]>0) {
    //
    //        BoundaryCylinderRepulsionExpenergy<<<blocksnthreadse[0], blocksnthreadse[1], (3 * blocksnthreadse[1]) * sizeof
    //                (floatingpoint), stream>>>
    //                          (coord, f, beadSet, krep, slen, nintvec, beListplane, params, gU_i, z,
    //                                  CUDAcommon::getCUDAvars().gculpritID, CUDAcommon::getCUDAvars().gculpritFF,
    //                                  CUDAcommon::getCUDAvars().gculpritinteraction, gFF, ginteraction);
    //        CUDAcommon::handleerror( cudaGetLastError() ,"BoundaryCylinderRepulsionExpenergy", "BoundaryCylinderRepulsionExp.cu");
    //    }
    
    if(blocksnthreadsez[1]>0) {
        std::cout<<"Total number of beads "<<Bead::getBeads().size()<<endl;
        auto boolvarvec = CUDAcommon::cudavars.backtrackbools;
        //        BoundaryCylinderRepulsionExpenergyz << < blocksnthreadsez[0], blocksnthreadsez[1], (6 * blocksnthreadsez[1]) *
        BoundaryCylinderRepulsionExpenergyz << < blocksnthreadsez[0], blocksnthreadsez[1], (0) *
        sizeof(floatingpoint), stream>> >(coord, f, beadSet, krep, slen, nintvec,
                                          beListplane, params, gU_i,
                                          CUDAcommon::cudavars.gpu_energyvec, z, CUDAcommon::getCUDAvars().gculpritID,
                                          CUDAcommon::getCUDAvars().gculpritFF,
                                          CUDAcommon::getCUDAvars()
                                          .gculpritinteraction, gFF,
                                          ginteraction, boolvarvec.at(0),
                                          boolvarvec.at(1));
        cudaDeviceSynchronize();
        
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
        CUDAcommon::handleerror(cudaGetLastError(),"BoundaryCylinderRepulsionExpenergyz", "BoundaryCylinderRepulsionExp.cu");
        return gU_sum;
    }
}
void BoundaryCylinderRepulsionExpIn::forces(floatingpoint *coord, floatingpoint *f, int *beadSet, floatingpoint *krep, floatingpoint *slen,
                                          int* nintvec, floatingpoint* beListplane, int *params) {
    if (blocksnthreadsf[1] > 0) {
        //        BoundaryCylinderRepulsionExpforces << < blocksnthreadsf[0], blocksnthreadsf[1], (3 * blocksnthreadsf[1]) *
        BoundaryCylinderRepulsionExpforces << < blocksnthreadsf[0], blocksnthreadsf[1], (0) *
        sizeof(floatingpoint), stream >> >(coord, f, beadSet, krep, slen, nintvec,
                                           beListplane, params);
        auto cvars = CUDAcommon::getCUDAvars();
        cvars.streamvec.push_back(&stream);
        CUDAcommon::cudavars = cvars;
        CUDAcommon::handleerror(cudaGetLastError(),"BoundaryCylinderRepulsionExpforces", "BoundaryCylinderRepulsionExp.cu");
    }
}
void BoundaryCylinderRepulsionExpIn::checkforculprit() {
    CUDAcommon::printculprit("BoundaryCylinderRepulsion","BoundaryCylinderRepulsionExp");
    Cylinder *c;
    BoundaryElement* be;
    be = BoundaryElement::getBoundaryElements()[CUDAcommon::getCUDAvars().culpritID[1]];
    cout << "Printing the culprit boundary element..." << endl;
    be->printSelf();
    c = Cylinder::getCylinders()[CUDAcommon::getCUDAvars().culpritID[0]];
    cout << "Printing the other culprit structure..." << endl;
    c->printSelf();
    exit(EXIT_FAILURE);
}
#endif
floatingpoint BoundaryCylinderRepulsionExpIn::energy(floatingpoint *coord, int *beadSet,
                                                   floatingpoint *krep, floatingpoint *slen, int *nneighbors) {
    
    int nb, nc;
    floatingpoint *coord1, R, r, U_i;
    floatingpoint U = 0.0;
    int Cumnc=0;
    auto beList = BoundaryElement::getBoundaryElements();
    nb = beList.size();
    
    //    for (int ib = 0; ib < nb; ib++) {
    //        auto be = beList[ib];
    //        be->printSelf();
    //    }
    
    for (int ib = 0; ib < nb; ib++) {
        
        auto be = beList[ib];
        nc = nneighbors[ib];
        
        for (int ic = 0; ic < nc; ic++) {
            
            coord1 = &coord[beadSet[Cumnc + ic]];
            r = be->distance(coord1);
            
            R = -r / slen[Cumnc + ic] + 100.0 / slen[Cumnc + ic];
            U_i = krep[Cumnc + ic] * exp(R);
            //            floatingpoint *var;
            //            var = new floatingpoint[4];
            //            be->elementeqn(var);
            //            std::cout<<"SL "<<Cumnc + ic<<" "<<krep[Cumnc + ic]<<" "<<R<<" Coord "<<coord1[0]<<" "<<coord1[1]<<" "
            //            <<coord1[2]<<" Plane "<<var[0] <<" "<<var[1]<<" "<<var[2]<<" "<<var[3]<<endl;
            //            delete var;
            //            std::cout<<r<<" "<<U_i<<endl;
            if (fabs(U_i) == numeric_limits<floatingpoint>::infinity()
                || U_i != U_i || U_i < -1.0) {
                
                //set culprit and return
                BoundaryInteractions::_boundaryElementCulprit = be;
                ///TODO
                //BoundaryInteractions::_otherCulprit;
                
                return -1;
            }
            U += U_i;
        }
        Cumnc += nc;
    }
    return U;
}

floatingpoint BoundaryCylinderRepulsionExpIn::energy(floatingpoint *coord, floatingpoint *f, int *beadSet,
                                                   floatingpoint *krep, floatingpoint *slen, int *nneighbors, floatingpoint d) {
    
    int nb, nc;
    floatingpoint *coord1, R, r, U_i;
    floatingpoint *force1;
    floatingpoint U = 0.0;
    int Cumnc=0;
    auto beList = BoundaryElement::getBoundaryElements();
    nb = beList.size();
    
    for (int ib = 0; ib < nb; ib++) {
        
        auto be = beList[ib];
        nc = nneighbors[ib];
        
        for(int ic = 0; ic < nc; ic++) {
            
            coord1 = &coord[beadSet[Cumnc + ic]];
            force1 = &f[beadSet[Cumnc + ic]];
            
            r = be->stretchedDistance(coord1, force1, d);
            
            R = -r / slen[Cumnc + ic] + 100.0 / slen[Cumnc + ic];
            
            //            std::cout<<r<<" "<<krep[Cumnc+ic]<<endl;
            U_i = krep[Cumnc + ic] * exp(R);
            //            floatingpoint *var;
            //            var = new floatingpoint[4];
            //            be->elementeqn(var);
            //            std::cout<<"SLZ "<<Cumnc + ic<<" "<<krep[Cumnc + ic]<<" "<<R<<" Coord "<<coord1[0]<<" "
            //                    <<coord1[1]<<" "<<coord1[2]<<" Plane "<<var[0] <<" "<<var[1]<<" "<<var[2]<<" "<<var[3]<<endl;
            //            delete var;
            if(fabs(U_i) == numeric_limits<floatingpoint>::infinity()
               || U_i != U_i || U_i < -1.0) {
                
                //set culprit and return
                BoundaryInteractions::_boundaryElementCulprit = be;
                ///TODO
                //BoundaryInteractions::_otherCulprit;
                
                return -1;
            }
            U += U_i;
        }
        Cumnc+=nc;
    }
    return U;
}



void BoundaryCylinderRepulsionExpIn::forces(floatingpoint *coord, floatingpoint *f, int *beadSet,
                                          floatingpoint *krep, floatingpoint *slen, int *nneighbors) {
    int nb, nc;
    floatingpoint *coord1, R, r, f0;
    floatingpoint *force1;
    //    floatingpoint *forcecopy;
    //    forcecopy = new floatingpoint[CGMethod::N];
    //    for(auto iter=0;iter<CGMethod::N;iter++)
    //        forcecopy[iter]=0.0;
    
    auto beList = BoundaryElement::getBoundaryElements();
    nb = beList.size();
    int Cumnc=0;
    
    for (int ib = 0; ib < nb; ib++) {
        
        auto be = beList[ib];
        nc = nneighbors[ib];
        for(int ic = 0; ic < nc; ic++) {
            coord1 = &coord[beadSet[ Cumnc + ic]];
            force1 = &f[beadSet[ Cumnc + ic]];
            r = be->distance(coord1);
            auto norm = be->normal(coord1);
            
            R = -r / slen[Cumnc + ic] + 100.0 / slen[Cumnc + ic];
            f0 = krep[Cumnc + ic] * exp(R)/ slen[Cumnc + ic];
            force1[0] += f0 *norm[0];
            force1[1] += f0 *norm[1];
            force1[2] += f0 *norm[2];
#ifdef CHECKFORCES_INF_NAN
            if(checkNaN_INF<floatingpoint>(force1, 0, 2)){
                cout<<"Boundary Cylinder Force becomes infinite. Printing data "<<endl;
                
                cout<<"Printing coords"<<endl;
                cout<<coord1[0]<<" "<<coord1[1]<<" "<<coord1[2]<<endl;
                cout<<"Printing force"<<endl;
                cout<<force1[0]<<" "<<force1[1]<<" "<<force1[2]<<endl;
                cout<<"Printing binary Coords"<<endl;
                printvariablebinary(coord1,0,2);
                cout<<"Printing binary Force"<<endl;
                printvariablebinary(force1,0,2);
                exit(EXIT_FAILURE);
            }
#endif
            
        }
        Cumnc+=nc;
    }}

floatingpoint BoundaryCylinderRepulsionExpIn::loadForces(floatingpoint r, floatingpoint kRep, floatingpoint screenLength) const {
    
    floatingpoint R = -r/screenLength + 100.0 / screenLength;
    return kRep * exp(R)/screenLength;
    
}
