
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
void CylinderExclVolume<CVolumeInteractionType>::vectorize(const FFCoordinateStartingIndex& si) {
    //count interactions
    nint = 0;

    for(auto ci : Cylinder::getCylinders()) {

        //do not calculate exvol for a non full length cylinder
#if defined(HYBRID_NLSTENCILLIST) || defined(SIMDBINDINGSEARCH)
        for (int ID = 0; ID < _HnlIDvec.size(); ID ++){
            auto neighbors = _HneighborList->getNeighborsstencil(_HnlIDvec[ID], ci);
            for(auto &cn : neighbors)
            {
                if(cn->getBranchingCylinder() == ci) continue;
                nint++;
            }
        }
#else
        auto neighbors = _neighborList->getNeighbors(ci);

        for(auto &cn : neighbors)
        {

            if(cn->getBranchingCylinder() == ci) continue;

            nint++;
        }
#endif
    }

    numInteractions = nint;
    CUDAcommon::tmin.numinteractions[8] += numInteractions;
//    std::cout<<"NINT1 "<<nint<<endl;
    beadSet = new int[n * nint];
    krep = new floatingpoint[nint];
    vecEqLength.resize(2 * nint);


    int nc = Cylinder::getCylinders().size();
    int i = 0;
    int Cumnc=0;
    for (i = 0; i < nc; i++) {
        auto ci = Cylinder::getCylinders()[i];

#if defined(HYBRID_NLSTENCILLIST) || defined(SIMDBINDINGSEARCH)
        
        for (int ID = 0; ID < _HnlIDvec.size(); ID ++){
            auto neighbors = _HneighborList->getNeighborsstencil(_HnlIDvec[ID], ci);
            int nn = neighbors.size();
            //        std::cout<<"Cylinder "<<i<<" "<<nn<<endl;
            for (int ni = 0; ni < nn; ni++) {
                
                auto cin = neighbors[ni];
                if(cin->getBranchingCylinder() == ci) continue;
                beadSet[n * (Cumnc)] = ci->getFirstBead()->getIndex()* 3 + si.bead;
                beadSet[n * (Cumnc) + 1] = ci->getSecondBead()->getIndex()* 3 + si.bead;
                beadSet[n * (Cumnc) + 2] = cin->getFirstBead()->getIndex()* 3 + si.bead;
                beadSet[n * (Cumnc) + 3] = cin->getSecondBead()->getIndex()* 3 + si.bead;

                vecEqLength[2 * Cumnc    ] = ci ->getMCylinder()->getEqLength();
                vecEqLength[2 * Cumnc + 1] = cin->getMCylinder()->getEqLength();
                
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
            }
        }
                
#else
        auto neighbors = _neighborList->getNeighbors(ci);
        
        int nn = neighbors.size();
//        std::cout<<"Cylinder "<<i<<" "<<nn<<endl;
        for (int ni = 0; ni < nn; ni++) {

            auto cin = neighbors[ni];
            if(cin->getBranchingCylinder() == ci) continue;
            beadSet[n * (Cumnc)] = ci->getFirstBead()->getIndex() * 3 + si.bead;
            beadSet[n * (Cumnc) + 1] = ci->getSecondBead()->getIndex() * 3 + si.bead;
            beadSet[n * (Cumnc) + 2] = cin->getFirstBead()->getIndex() * 3 + si.bead;
            beadSet[n * (Cumnc) + 3] = cin->getSecondBead()->getIndex() * 3 + si.bead;

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
#endif
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
//   if(blocksnthreads[0]==1) blocksnthreads.push_back( 32*(int(numInteractions/32 +1)) );
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
    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_krep, numInteractions * sizeof(floatingpoint)),"cuda data transfer",
                            "CylinderExclVolume.cu");
    CUDAcommon::handleerror(cudaMemcpyAsync(gpu_krep, krep, numInteractions * sizeof
                            (floatingpoint), cudaMemcpyHostToDevice, stream),
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
    chrono::duration<floatingpoint> elapsed_run(tend - tbegin);
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
floatingpoint CylinderExclVolume<CVolumeInteractionType>::computeEnergy(floatingpoint *coord) {


    floatingpoint U_ii=0.0f;

#ifdef CUDATIMETRACK
    chrono::high_resolution_clock::time_point tbegin, tend;
#endif
#ifdef CUDAACCL

#ifdef CUDATIMETRACK
    tbegin = chrono::high_resolution_clock::now();
#endif

    //has to be changed to accomodate aux force
    floatingpoint * gpu_coord=CUDAcommon::getCUDAvars().gpu_coord;
    floatingpoint * gpu_force = CUDAcommon::getCUDAvars().gpu_force;
    floatingpoint * gpu_d = CUDAcommon::getCUDAvars().gpu_lambda;

//    if(d == 0.0){
//        gU_i=_FFType.energy(gpu_coord, gpu_force, gpu_beadSet, gpu_krep, gpu_params);
//    }
//    else{
        gU_i=_FFType.energy(gpu_coord, gpu_force, gpu_beadSet, gpu_krep, gpu_d, gpu_params);
//    }

#ifdef CUDATIMETRACK
//    CUDAcommon::handleerror(cudaDeviceSynchronize(),"CylinderExclVolume.cu", "computeEnergy");
    tend= chrono::high_resolution_clock::now();
    chrono::duration<floatingpoint> elapsed_run(tend - tbegin);
    CUDAcommon::cudatime.TveccomputeE.push_back(elapsed_run.count());
    CUDAcommon::cudatime.TcomputeE += elapsed_run.count();
    CUDAcommon::cudatime.TcomputeEiter += elapsed_run.count();
#endif
#endif
#ifdef SERIAL
#ifdef CUDATIMETRACK
    tbegin = chrono::high_resolution_clock::now();
#endif

    U_ii = _FFType.energy(coord, beadSet, krep, vecEqLength.data(), numInteractions);

#ifdef CUDATIMETRACK
    floatingpoint U_i[1];
    floatingpoint *gU_i;
    tend= chrono::high_resolution_clock::now();
    chrono::duration<floatingpoint> elapsed_runs(tend - tbegin);
    CUDAcommon::serltime.TveccomputeE.push_back(elapsed_runs.count());
    CUDAcommon::serltime.TcomputeE += elapsed_runs.count();
    CUDAcommon::serltime.TcomputeEiter += elapsed_runs.count();
#endif

#endif
    return U_ii;
}

template <class CVolumeInteractionType>
void CylinderExclVolume<CVolumeInteractionType>::computeForces(floatingpoint *coord, floatingpoint *f) {
#ifdef CUDATIMETRACK
    chrono::high_resolution_clock::time_point tbegin, tend;
    tbegin = chrono::high_resolution_clock::now();
#endif
#ifdef CUDAACCL
    //has to be changed to accomodate aux force
    floatingpoint * gpu_coord=CUDAcommon::getCUDAvars().gpu_coord;
    floatingpoint * gpu_force;
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
    chrono::duration<floatingpoint> elapsed_run(tend - tbegin);
    CUDAcommon::cudatime.TveccomputeF.push_back(elapsed_run.count());
    CUDAcommon::cudatime.TcomputeF += elapsed_run.count();
    tbegin = chrono::high_resolution_clock::now();
#endif
#ifdef SERIAL
    _FFType.forces(coord, f, beadSet, krep, vecEqLength.data(), numInteractions);
#endif
#ifdef CUDATIMETRACK
    tend= chrono::high_resolution_clock::now();
    chrono::duration<floatingpoint> elapsed_runs(tend - tbegin);
    CUDAcommon::serltime.TveccomputeF.push_back(elapsed_runs.count());
    CUDAcommon::serltime.TcomputeF += elapsed_runs.count();
#endif
#ifdef DETAILEDOUTPUT
    floatingpoint maxF = 0.0;
    floatingpoint mag = 0.0;
    for(int i = 0; i < CGMethod::N/3; i++) {
        mag = 0.0f;
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
template floatingpoint CylinderExclVolume<CylinderExclVolRepulsion>::computeEnergy(floatingpoint *coord);
template void CylinderExclVolume<CylinderExclVolRepulsion>::computeForces(floatingpoint *coord, floatingpoint *f);
template void CylinderExclVolume<CylinderExclVolRepulsion>::vectorize(const FFCoordinateStartingIndex&);
template void CylinderExclVolume<CylinderExclVolRepulsion>::deallocate();
