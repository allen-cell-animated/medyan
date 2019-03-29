
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

#include "CylinderExclVolume.h"

#include "CylinderExclVolRepulsion.h"

#include "Cylinder.h"
#include "Bead.h"

#include "MathFunctions.h"
#include "cross_check.h"
#include "CGMethod.h"
#ifdef CUDAACCL
#include "nvToolsExt.h"
#endif


using namespace mathfunc;

template <class CVolumeInteractionType>
int CylinderExclVolume<CVolumeInteractionType>::numInteractions;

template <class CVolumeInteractionType>
void CylinderExclVolume<CVolumeInteractionType>::vectorize() {
    //count interactions
    nint = 0;

    for(auto ci : Cylinder::getCylinders()) {
        if(!ci->isFullLength()) continue;
        //do not calculate exvol for a non full length cylinder
#ifdef HYBRID_NLSTENCILLIST
        auto neighbors = _HneighborList->getNeighborsstencil(_HnlID, ci);
#else
        auto neighbors = _neighborList->getNeighbors(ci);
#endif
        for(auto &cn : neighbors)
        {
            if(!cn->isFullLength()||
               cn->getBranchingCylinder() == ci) continue;
            nint++;
        }
    }

    numInteractions = nint;
//    std::cout<<"NINT1 "<<nint<<endl;
    beadSet = new int[n * nint];
    krep = new double[nint];


    int nc = Cylinder::getCylinders().size();
    int i = 0;
    int Cumnc=0;
    for (i = 0; i < nc; i++) {
        auto ci = Cylinder::getCylinders()[i];
        if(!ci->isFullLength()) continue;
#ifdef HYBRID_NLSTENCILLIST
        auto neighbors = _HneighborList->getNeighborsstencil(_HnlID, ci);
#else
        auto neighbors = _neighborList->getNeighbors(ci);
#endif
        int nn = neighbors.size();
//        std::cout<<"Cylinder "<<i<<" "<<nn<<endl;
        for (int ni = 0; ni < nn; ni++) {

            auto cin = neighbors[ni];
            if(!cin->isFullLength()||
               cin->getBranchingCylinder() == ci) continue;
            beadSet[n * (Cumnc)] = ci->getFirstBead()->_dbIndex;
            beadSet[n * (Cumnc) + 1] = ci->getSecondBead()->_dbIndex;
            beadSet[n * (Cumnc) + 2] = cin->getFirstBead()->_dbIndex;
            beadSet[n * (Cumnc) + 3] = cin->getSecondBead()->_dbIndex;
            
            //Get KRepuls based on filament type
            if(ci->getType() != cin->getType()){
                auto ki = ci->getMCylinder()->getExVolConst();
                auto kin = cin->getMCylinder()->getExVolConst();
                krep[Cumnc] = max(ki, kin);
            }
            else{
                krep[Cumnc] = ci->getMCylinder()->getExVolConst();
            }
            
            Cumnc++;
            //std::cout<<"CV"<<ci->getID()<<" "<<cin->getID()<<endl;
        }
    }
    //CUDA
#ifdef CUDAACCL
#ifdef CUDATIMETRACK
    chrono::high_resolution_clock::time_point tbegin, tend;
    tbegin = chrono::high_resolution_clock::now();
#endif

//    cudaEvent_t start, stop;
//    CUDAcommon::handleerror(cudaEventCreate( &start));
//    CUDAcommon::handleerror(cudaEventCreate( &stop));
//    CUDAcommon::handleerror(cudaEventRecord( start, 0));

//    blocksnthreads.push_back(int(numInteractions/THREADSPERBLOCK + 1));
//    if(blocksnthreads[0]==1) blocksnthreads.push_back( numInteractions);
////    if(blocksnthreads[0]==1) blocksnthreads.push_back( 32*(int(numInteractions/32 +1)) );
//    else blocksnthreads.push_back(THREADSPERBLOCK);

    //CUDA stream create
    if(stream == NULL || !(CUDAcommon::getCUDAvars().conservestreams))
        CUDAcommon::handleerror(cudaStreamCreate(&stream),"cuda stream", "CylinderExclVolume.cu");
    _FFType.optimalblocksnthreads(numInteractions, stream);

    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_beadSet, n * numInteractions * sizeof(int)),"cuda data transfer",
                            "CylinderExclVolume.cu");
    CUDAcommon::handleerror(cudaMemcpyAsync(gpu_beadSet, beadSet, n * numInteractions *
                                    sizeof(int), cudaMemcpyHostToDevice, stream),
                            "cuda data transfer", "CylinderExclVolume.cu");
    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_krep, numInteractions * sizeof(double)),"cuda data transfer",
                            "CylinderExclVolume.cu");
    CUDAcommon::handleerror(cudaMemcpyAsync(gpu_krep, krep, numInteractions * sizeof
                            (double), cudaMemcpyHostToDevice, stream),
                            "cuda data transfer", "CylinderExclVolume.cu");
    vector<int> params;
    params.push_back(int(n));
    params.push_back(numInteractions);
    params.push_back(CUDAcommon::cudavars.offset_E);
//    std::cout<<"CUDA exvol offsetE "<<CUDAcommon::cudavars.offset_E<<endl;
    //set offset
    CUDAcommon::cudavars.offset_E += nint;
//    std::cout<<"offset "<<getName()<<" "<<CUDAcommon::cudavars.offset_E<<endl;
    //TODO make sure not using cudafree here is fine.
    if(gpu_params != NULL )
        CUDAcommon::handleerror(cudaFree(gpu_params),"cudaFree", "CylinderExclVolume.cu");
    if(nint > 0) {
        CUDAcommon::handleerror(cudaMalloc((void **) &gpu_params, 3 * sizeof(int)), "cuda"
                                        " data transfer",
                                "CylinderExclVolume.cu");
        CUDAcommon::handleerror(cudaMemcpyAsync(gpu_params, params.data(), 3 * sizeof(int),
                                           cudaMemcpyHostToDevice, stream),
                                "cuda data transfer", "CylinderExclVolume.cu");
    }
#ifdef CUDATIMETRACK
//    CUDAcommon::handleerror(cudaDeviceSynchronize(),"CylinderExclVolume.cu",
//                            "vectorizeFF");
    tend= chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed_run(tend - tbegin);
    CUDAcommon::cudatime.TvecvectorizeFF.push_back(elapsed_run.count());
    CUDAcommon::cudatime.TvectorizeFF += elapsed_run.count();
#endif
#endif
    //
}


template <class CVolumeInteractionType>
void CylinderExclVolume<CVolumeInteractionType>::deallocate() {

    delete [] beadSet;
    delete [] krep;
#ifdef CUDAACCL
    if(nint > 0) {
        if(!(CUDAcommon::getCUDAvars().conservestreams))
            CUDAcommon::handleerror(cudaStreamDestroy(stream),"cuda stream", "CylinderExclVolume.cu");
        _FFType.deallocate();
        CUDAcommon::handleerror(cudaFree(gpu_beadSet), "cudaFree", "CylinderExclVolume.cu");
        CUDAcommon::handleerror(cudaFree(gpu_krep), "cudaFree", "CylinderExclVolume.cu");
        CUDAcommon::handleerror(cudaFree(gpu_params), "cudaFree", "CylinderExclVolume.cu");
        gpu_beadSet = NULL;
        gpu_krep = NULL;
        gpu_params = NULL;
    }
#endif
}


template <class CVolumeInteractionType>
double CylinderExclVolume<CVolumeInteractionType>::computeEnergy(double *coord, double *f, double d) {

    double U_i[1];
    double U_ii=0.0;
    double *gU_i;
#ifdef CUDATIMETRACK
    chrono::high_resolution_clock::time_point tbegin, tend;
#endif
#ifdef CUDAACCL

#ifdef CUDATIMETRACK
    tbegin = chrono::high_resolution_clock::now();
#endif

    //has to be changed to accomodate aux force
    double * gpu_coord=CUDAcommon::getCUDAvars().gpu_coord;
    double * gpu_force = CUDAcommon::getCUDAvars().gpu_force;
    double * gpu_d = CUDAcommon::getCUDAvars().gpu_lambda;

//    if(d == 0.0){
//        gU_i=_FFType.energy(gpu_coord, gpu_force, gpu_beadSet, gpu_krep, gpu_params);
//    }
//    else{
        gU_i=_FFType.energy(gpu_coord, gpu_force, gpu_beadSet, gpu_krep, gpu_d, gpu_params);
//    }

#ifdef CUDATIMETRACK
//    CUDAcommon::handleerror(cudaDeviceSynchronize(),"CylinderExclVolume.cu", "computeEnergy");
    tend= chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed_run(tend - tbegin);
    CUDAcommon::cudatime.TveccomputeE.push_back(elapsed_run.count());
    CUDAcommon::cudatime.TcomputeE += elapsed_run.count();
    CUDAcommon::cudatime.TcomputeEiter += elapsed_run.count();
#endif
#endif
#ifdef SERIAL
#ifdef CUDATIMETRACK
    tbegin = chrono::high_resolution_clock::now();
#endif
    if (d == 0.0)
        U_ii = _FFType.energy(coord, f, beadSet, krep);
    else
        U_ii = _FFType.energy(coord, f, beadSet, krep, d);

#ifdef CUDATIMETRACK
    tend= chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed_runs(tend - tbegin);
    CUDAcommon::serltime.TveccomputeE.push_back(elapsed_runs.count());
    CUDAcommon::serltime.TcomputeE += elapsed_runs.count();
    CUDAcommon::serltime.TcomputeEiter += elapsed_runs.count();
#endif

#endif
    return U_ii;
}

template <class CVolumeInteractionType>
void CylinderExclVolume<CVolumeInteractionType>::computeForces(double *coord, double *f) {
#ifdef CUDATIMETRACK
    chrono::high_resolution_clock::time_point tbegin, tend;
    tbegin = chrono::high_resolution_clock::now();
#endif
#ifdef CUDAACCL
    //has to be changed to accomodate aux force
    double * gpu_coord=CUDAcommon::getCUDAvars().gpu_coord;
    double * gpu_force;
    if(cross_checkclass::Aux) {

        gpu_force = CUDAcommon::getCUDAvars().gpu_forceAux;
        _FFType.forces(gpu_coord, gpu_force, gpu_beadSet, gpu_krep, gpu_params);

    }
    else {

        gpu_force = CUDAcommon::getCUDAvars().gpu_force;
        _FFType.forces(gpu_coord, gpu_force, gpu_beadSet, gpu_krep, gpu_params);
    }
#endif
#ifdef CUDATIMETRACK
    tend= chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed_run(tend - tbegin);
    CUDAcommon::cudatime.TveccomputeF.push_back(elapsed_run.count());
    CUDAcommon::cudatime.TcomputeF += elapsed_run.count();
    tbegin = chrono::high_resolution_clock::now();
#endif
#ifdef SERIAL
    _FFType.forces(coord, f, beadSet, krep);
#endif
#ifdef CUDATIMETRACK
    tend= chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed_runs(tend - tbegin);
    CUDAcommon::serltime.TveccomputeF.push_back(elapsed_runs.count());
    CUDAcommon::serltime.TcomputeF += elapsed_runs.count();
#endif
#ifdef DETAILEDOUTPUT
    double maxF = 0.0;
    double mag = 0.0;
    for(int i = 0; i < CGMethod::N/3; i++) {
        mag = 0.0;
        for(int j = 0; j < 3; j++)
            mag += f[3 * i + j]*f[3 * i + j];
        mag = sqrt(mag);
//        std::cout<<"SL "<<i<<" "<<mag*mag<<" "<<forceAux[3 * i]<<" "<<forceAux[3 * i + 1]<<" "<<forceAux[3 * i +
//                2]<<endl;
        if(mag > maxF) maxF = mag;
    }
    std::cout<<"max "<<getName()<<" "<<maxF<<endl;
#endif
}

///Template specializations
template double CylinderExclVolume<CylinderExclVolRepulsion>::computeEnergy(double *coord, double *f, double d);
template void CylinderExclVolume<CylinderExclVolRepulsion>::computeForces(double *coord, double *f);
template void CylinderExclVolume<CylinderExclVolRepulsion>::vectorize();
template void CylinderExclVolume<CylinderExclVolRepulsion>::deallocate();


