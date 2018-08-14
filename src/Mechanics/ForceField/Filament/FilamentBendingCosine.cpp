
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

#include "FilamentBendingCosine.h"
#include "FilamentBending.h"
#include "FilamentBendingCosineCUDA.h"

#include "Filament.h"
#include "Cylinder.h"
#include "Bead.h"

#include "MathFunctions.h"

#ifdef CUDAACCL
#include <cuda.h>
#include <cuda_runtime.h>
#include "nvToolsExt.h"
#endif

using namespace mathfunc;
#ifdef CUDAACCL
void FilamentBendingCosine::deallocate(){
    if(!(CUDAcommon::getCUDAvars().conservestreams))
        CUDAcommon::handleerror(cudaStreamDestroy(stream));
    CUDAcommon::handleerror(cudaFree(gU_i));
    CUDAcommon::handleerror(cudaFree(gU_sum));
    CUDAcommon::handleerror(cudaFree(gFF));
    CUDAcommon::handleerror(cudaFree(ginteraction));
}
void FilamentBendingCosine::optimalblocksnthreads( int nint){
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
                                                       FilamentBendingCosineenergy, blockToSmemFB, 0);
        blocksnthreadse.push_back((nint + blockSize - 1) / blockSize);
        blocksnthreadse.push_back(blockSize);
        blockSize = 0;

        cudaOccupancyMaxPotentialBlockSizeVariableSMem(&minGridSize, &blockSize,
                                                       FilamentBendingCosineenergyz, blockToSmemFB2, 0);
        blocksnthreadsez.push_back((nint + blockSize - 1) / blockSize);
        blocksnthreadsez.push_back(blockSize);
        blockSize = 0;

        cudaOccupancyMaxPotentialBlockSizeVariableSMem(&minGridSize, &blockSize,
                                                       FilamentBendingCosineforces, blockToSmemFB, 0);
        blocksnthreadsf.push_back((nint + blockSize - 1) / blockSize);
        blocksnthreadsf.push_back(blockSize);
        //get addition vars
        bntaddvec2.clear();
        bntaddvec2 = getaddred2bnt(nint);
        CUDAcommon::handleerror(cudaMalloc((void **) &gU_i, bntaddvec2.at(0)*sizeof(double)));
        CUDAcommon::handleerror(cudaMemset(gU_i, 0, bntaddvec2.at(0) * sizeof(double)));
//        CUDAcommon::handleerror(cudaMalloc((void **) &gU_i, nint*sizeof(double)));
        CUDAcommon::handleerror(cudaMalloc((void **) &gU_sum, sizeof(double)));

        char a[] = "FilamentFF";
        char b[] = "Filament Bending Cosine";
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
double* FilamentBendingCosine::energy(double *coord, double *f, int *beadSet,
                                      double *kbend, double *eqt, int *params) {
//    if(blocksnthreadse[1]>0) {
//        FilamentBendingCosineenergy<<<blocksnthreadse[0], blocksnthreadse[1], (9 * blocksnthreadse[1]) * sizeof
//                (double), stream>>> (coord, f, beadSet, kbend, eqt, params, gU_i, CUDAcommon::getCUDAvars()
//                .gculpritID,
//                CUDAcommon::getCUDAvars().gculpritFF,
//                CUDAcommon::getCUDAvars().gculpritinteraction, gFF, ginteraction);
//        auto cvars = CUDAcommon::getCUDAvars();
//        cvars.streamvec.push_back(&stream);
//        CUDAcommon::cudavars = cvars;
//        CUDAcommon::handleerror( cudaGetLastError() ,"FilamentBendingCosineenergy", "FilamentBendingCosine.cu");
//        double* gpu_Utot = CUDAcommon::getCUDAvars().gpu_energy;
//        addvector<<<1,1,0,stream>>>(gU_i,params, gU_sum, gpu_Utot);
//        CUDAcommon::handleerror( cudaGetLastError() ,"FilamentBendingCosineenergy", "FilamentBendingCosine.cu");
//        return gU_sum;}
//    else
//        return NULL;
}


double* FilamentBendingCosine::energy(double *coord, double *f, int *beadSet,
                                      double *kbend, double *eqt, double *z, int *params) {
//    nvtxRangePushA("E_wait");
//    CUDAcommon::handleerror(cudaStreamWaitEvent(stream, *(CUDAcommon::getCUDAvars().event), 0));
//    nvtxRangePop();
    if(blocksnthreadse[1]>0) {
        FilamentBendingCosineenergy<<<blocksnthreadse[0], blocksnthreadse[1], (9 * blocksnthreadse[1]) * sizeof
                (double), stream>>> (coord, f, beadSet, kbend, eqt, params, gU_i, z, CUDAcommon::getCUDAvars()
                .gculpritID,
                CUDAcommon::getCUDAvars().gculpritFF,
                CUDAcommon::getCUDAvars().gculpritinteraction, gFF, ginteraction);
        CUDAcommon::handleerror(cudaGetLastError(),"FilamentBendingCosineenergy", "FilamentBendingCosine.cu");
//        auto cvars = CUDAcommon::getCUDAvars();
//        cvars.streamvec.push_back(&stream);
//        CUDAcommon::cudavars = cvars;
//        CUDAcommon::handleerror( cudaGetLastError() ,"FilamentBendingCosineenergy", "FilamentBendingCosine.cu");
//        double* gpu_Utot = CUDAcommon::getCUDAvars().gpu_energy;
//        addvector<<<1,1,0,stream>>>(gU_i,params, gU_sum, gpu_Utot);
//        CUDAcommon::handleerror( cudaGetLastError() ,"FilamentBendingCosineenergy", "FilamentBendingCosine.cu");
//        return gU_sum;
    }

    if(blocksnthreadsez[1]>0) {
        FilamentBendingCosineenergyz << < blocksnthreadsez[0], blocksnthreadsez[1], (18 * blocksnthreadsez[1]) *
                                          sizeof(double), stream>> > (coord, f, beadSet, kbend, eqt, params, gU_i, z,
                                          CUDAcommon::getCUDAvars().gculpritID,
                                          CUDAcommon::getCUDAvars().gculpritFF,
                                          CUDAcommon::getCUDAvars().gculpritinteraction, gFF, ginteraction);
        CUDAcommon::handleerror(cudaGetLastError(),"FilamentBendingCosineenergyz", "FilamentBendingCosine.cu");
//        auto cvars = CUDAcommon::getCUDAvars();
//        cvars.streamvec.push_back(&stream);
//        CUDAcommon::cudavars = cvars;
//        CUDAcommon::handleerror(cudaGetLastError(),"FilamentBendingCosineenergyz", "FilamentBendingCosine.cu");
//        double* gpu_Utot = CUDAcommon::getCUDAvars().gpu_energy;
//
//        addvector<<<1,1,0,stream>>>(gU_i,params, gU_sum, gpu_Utot);
//        CUDAcommon::handleerror(cudaGetLastError(),"FilamentBendingCosineenergyz", "FilamentBendingCosine.cu");
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

void FilamentBendingCosine::forces(double *coord, double *f, int *beadSet,
                                   double *kbend, double *eqt, int *params){
    if(blocksnthreadsf[1]>0) {
        FilamentBendingCosineforces << < blocksnthreadsf[0], blocksnthreadsf[1], (9 * blocksnthreadsf[1]) *
                                                                                 sizeof(double), stream >> > (coord, f, beadSet, kbend, eqt, params);
        auto cvars = CUDAcommon::getCUDAvars();
        cvars.streamvec.push_back(&stream);
        CUDAcommon::cudavars = cvars;
        CUDAcommon::handleerror(cudaGetLastError(),"FilamentBendingCosineforces", "FilamentBendingCosine.cu");
//        CUDAcommon::handleerror(cudaDeviceSynchronize());
//        double U_i[Bead::getBeads().size() - 2 * Filament::getFilaments().size()];
//        CUDAcommon::handleerror(cudaMemcpy(U_i, eqt, (Bead::getBeads().size() - 2 * Filament::getFilaments().size()) *
//                                                                               sizeof(double), cudaMemcpyDeviceToHost));
//        for(auto i = 0; i< Bead::getBeads().size() - 2 * Filament::getFilaments().size();i ++){
//            std::cout<<U_i[i]<<endl;
//        }
//        std::cout<<endl;
    }
}

void FilamentBendingCosine::checkforculprit() {
    CUDAcommon::printculprit("FilamentBending","FilamentBendingCosine");
    Filament* fil;
    int i = 0;
    bool found = false;
    for (auto f: Filament::getFilaments()) {

        if (f->getCylinderVector().size() > 1){
            i = i + 2 * f->getCylinderVector().size() - 2;
            if(i > CUDAcommon::getCUDAvars().culpritID[0] && !found){
                found = true;
                fil = (Filament*)(Cylinder::getCylinders()[i]->getParent());
            }
        }
    }
    cout<<"Printing culprit Filament information."<<endl;
    fil->printSelf();
    exit(EXIT_FAILURE);
}
#endif
double FilamentBendingCosine::energy(double *coord, double *f, int *beadSet,
                                     double *kbend, double *eqt){

    int n = FilamentBending<FilamentBendingCosine>::n;
    int nint = (Bead::getBeads().size() - 2 * Filament::getFilaments().size());

    double *coord1, *coord2, *coord3, U_i, L1, L2, L1L2, l1l2, phi, dPhi;

    double U = 0.0;

    for(int i = 0; i < nint; i += 1) {

        coord1 = &coord[3 * beadSet[n * i]];
        coord2 = &coord[3 * beadSet[n * i + 1]];
        coord3 = &coord[3 * beadSet[n * i + 2]];

        L1 = sqrt(scalarProduct(coord1, coord2,
                                coord1, coord2));
        L2 = sqrt(scalarProduct(coord2, coord3,
                                coord2, coord3));

        L1L2 = L1*L2;
        l1l2 = scalarProduct(coord1, coord2,
                             coord2, coord3);

        phi = safeacos(l1l2 / L1L2);
        dPhi = phi-eqt[i];

        U_i = kbend[i] * ( 1 - cos(dPhi) );

        if(fabs(U_i) == numeric_limits<double>::infinity()
           || U_i != U_i || U_i < -1.0) {

            //set culprit and return TODO
            //FilamentInteractions::_filamentCulprit = Filament::getFilaments()[i];

            return -1;
        }

        U += U_i;
    }

    return U;
}

double FilamentBendingCosine::energy(double *coord, double *f, int *beadSet,
                                     double *kbend, double *eqt, double d ){

    int n = FilamentBending<FilamentBendingCosine>::n;
    int nint =  (Bead::getBeads().size() - 2 * Filament::getFilaments().size());

    double *coord1, *coord2, *coord3, *force1, *force2, *force3, U_i, L1, L2, L1L2, l1l2, phi, dPhi;

    double U = 0.0;

    for(int i = 0; i < nint; i += 1) {

        coord1 = &coord[3 * beadSet[n * i]];
        coord2 = &coord[3 * beadSet[n * i + 1]];
        coord3 = &coord[3 * beadSet[n * i + 2]];

        force1 = &f[3 * beadSet[n * i]];
        force2 = &f[3 * beadSet[n * i + 1]];
        force3 = &f[3 * beadSet[n * i + 2]];


        L1 = sqrt(scalarProductStretched(coord1, force1, coord2, force2,
                                         coord1, force1, coord2, force2, d));
        L2 = sqrt(scalarProductStretched(coord2, force2, coord3, force3,
                                         coord2, force2, coord3, force3, d));

        L1L2 = L1*L2;
        l1l2 = scalarProductStretched(coord1, force1, coord2, force2,
                                      coord2, force2, coord3, force3, d);

        phi = safeacos(l1l2 / L1L2);
        dPhi = phi-eqt[i];

        U_i = kbend[i] * ( 1 - cos(dPhi) );

        if(fabs(U_i) == numeric_limits<double>::infinity()
           || U_i != U_i || U_i < -1.0) {

            //set culprit and return TODO
            //FilamentInteractions::_filamentCulprit = Filament::getFilaments()[i];

            return -1;
        }

        U += U_i;
    }

    return U;
}

void FilamentBendingCosine::forces(double *coord, double *f, int *beadSet,
                                   double *kbend, double *eqt){

    int n = FilamentBending<FilamentBendingCosine>::n;
    int nint =  (Bead::getBeads().size() - 2 * Filament::getFilaments().size());

    double *coord1, *coord2, *coord3, *force1, *force2, *force3,
            L1, L2, l1l2, invL1, invL2, A,B,C, phi, dPhi, k;

    for(int i = 0; i < nint; i += 1) {

        coord1 = &coord[3 * beadSet[n * i]];
        coord2 = &coord[3 * beadSet[n * i + 1]];
        coord3 = &coord[3 * beadSet[n * i + 2]];

        force1 = &f[3 * beadSet[n * i]];
        force2 = &f[3 * beadSet[n * i + 1]];
        force3 = &f[3 * beadSet[n * i + 2]];

        L1 = sqrt(scalarProduct(coord1, coord2,
                                coord1, coord2));
        L2 = sqrt(scalarProduct(coord2, coord3,
                                coord2, coord3));

        l1l2 = scalarProduct(coord1, coord2,
                             coord2, coord3);

        invL1 = 1/L1;
        invL2 = 1/L2;
        A = invL1*invL2;
        B = l1l2*invL1*A*A*L2;
        C = l1l2*invL2*A*A*L1;

        if (areEqual(eqt[i], 0.0)) k = kbend[i];

        else{
            phi = safeacos(l1l2 *A);
            dPhi = phi-eqt[i];

            k = kbend[i] * sin(dPhi)/sin(phi);
            cout<< "Watch out! Normally we set equalibrium theta to 0.0" << endl;
        }
        //force on i-1, f = k*(-A*l2 + B*l1):
        force1[0] +=  k * ((-coord3[0] + coord2[0])*A +
                           (coord2[0] - coord1[0])*B );
        force1[1] +=  k * ((-coord3[1] + coord2[1])*A +
                           (coord2[1] - coord1[1])*B );
        force1[2] +=  k * ((-coord3[2] + coord2[2])*A +
                           (coord2[2] - coord1[2])*B );


        //force on i, f = k*(A*(l1-l2) - B*l1 + C*l2):
        force2[0] +=  k *( (coord3[0] - 2*coord2[0] + coord1[0])*A -
                           (coord2[0] - coord1[0])*B +
                           (coord3[0] - coord2[0])*C );

        force2[1] +=  k *( (coord3[1] - 2*coord2[1] + coord1[1])*A -
                           (coord2[1] - coord1[1])*B +
                           (coord3[1] - coord2[1])*C );

        force2[2] +=  k *( (coord3[2] - 2*coord2[2] + coord1[2])*A -
                           (coord2[2] - coord1[2])*B +
                           (coord3[2] - coord2[2])*C );

        //force on i-1, f = k*(A*l - B*l2):
        force3[0] +=  k *( (coord2[0] - coord1[0])*A -
                           (coord3[0] - coord2[0])*C );

        force3[1] +=  k *( (coord2[1] - coord1[1])*A -
                           (coord3[1] - coord2[1])*C );

        force3[2] +=  k *( (coord2[2] - coord1[2])*A -
                           (coord3[2] - coord2[2])*C );
        
        double f1sq = force1[0] * force1[0] + force1[1] * force1[1] + force1[2] * force1[2];
        double f2sq = force2[0] * force2[0] + force2[1] * force2[1] + force2[2] * force2[2];
        double f3sq = force3[0] * force3[0] + force3[1] * force3[1] + force3[2] * force3[2];
        
        if(f1sq > 1e8 || f2sq > 1e8 || f3sq > 1e8){
            cout<<"High bending cosine force!" << endl;
            cout<<"coord1 = " << coord1[0] << " " << coord1[1] << " " << coord1[2] <<endl;
            cout<< "force1 = " << sqrt(f1sq) << endl;
            cout<<"coord2 = " << coord2[0] << " " << coord2[1] << " " << coord2[2] <<endl;
            cout<< "force2 = " << sqrt(f2sq) << endl;
            cout<<"coord3 = " << coord3[0] << " " << coord3[1] << " " << coord3[2] <<endl;
            cout<< "force3 = " << sqrt(f3sq) << endl;
            cout << "L1 = " << L1 << ", L2 = " << L2 << ", l1l2 = " << l1l2 <<endl;
            cout << "A = " << A << ", B = " << B << ", C = " << C <<endl;
            cout << "k = " << k  << ", eqt = " << eqt[i] << endl;
            cout << "nint = " << nint<<endl;
            for (auto b: Bead::getBeads()){
                if(b->_dbIndex == beadSet[beadSet[n * i]]){
                    //Then access all elements of beads that you want to access
                    b->printSelf();
                }
                if(b->_dbIndex == beadSet[beadSet[n * i] +1]){
                    //Then access all elements of beads that you want to access
                    b->printSelf();
                }
                if(b->_dbIndex == beadSet[beadSet[n * i] +2]){
                    //Then access all elements of beads that you want to access
                    b->printSelf();
                }
            }
        }
        
//        double f1[3], f2[3], f3[3];
//        f1[0] =  k * ((-coord3[0] + coord2[0])*A +
//                           (coord2[0] - coord1[0])*B );
//        f1[1] =  k * ((-coord3[1] + coord2[1])*A +
//                           (coord2[1] - coord1[1])*B );
//        f1[2] =  k * ((-coord3[2] + coord2[2])*A +
//                           (coord2[2] - coord1[2])*B );
//
//
//        //force on i, f = k*(A*(l1-l2) - B*l1 + C*l2):
//        f2[0] =  k *( (coord3[0] - 2*coord2[0] + coord1[0])*A -
//                           (coord2[0] - coord1[0])*B +
//                           (coord3[0] - coord2[0])*C );
//
//        f2[1] =  k *( (coord3[1] - 2*coord2[1] + coord1[1])*A -
//                           (coord2[1] - coord1[1])*B +
//                           (coord3[1] - coord2[1])*C );
//
//        f2[2] =  k *( (coord3[2] - 2*coord2[2] + coord1[2])*A -
//                           (coord2[2] - coord1[2])*B +
//                           (coord3[2] - coord2[2])*C );
//
//        //force on i-1, f = k*(A*l - B*l2):
//        f3[0] =  k *( (coord2[0] - coord1[0])*A -
//                           (coord3[0] - coord2[0])*C );
//
//        f3[1] =  k *( (coord2[1] - coord1[1])*A -
//                           (coord3[1] - coord2[1])*C );
//
//        f3[2] =  k *( (coord2[2] - coord1[2])*A -
//                           (coord3[2] - coord2[2])*C );
//        std::cout<<i<<" "<<f1[0]<<" "<<f1[1]<<" "<<f1[2]<<" "<<f2[0]<<" "<<f2[1]<<" "<<f2[2]<<" "<<f3[0]<<" "
//                ""<<f3[1]<<" "<<f3[2]<<endl;

//               std::cout<<"BENDING "<<force1[0]<<" "<<force1[1]<<" "<<force1[2]<<" "<<force2[0]<<" "<<force2[1]<<" "<<force2[2]<<" "<<force3[0]<<" "<<force3[1]<<" "<<force3[2]<<endl;
    }
}
