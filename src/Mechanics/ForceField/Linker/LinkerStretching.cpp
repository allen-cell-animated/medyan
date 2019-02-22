
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

#include "LinkerStretching.h"

#include "LinkerStretchingHarmonic.h"

#include "Cylinder.h"
#include "Linker.h"
#include "Bead.h"
#include "cross_check.h"
#include "CGMethod.h"
#ifdef CUDAACCL
#include "nvToolsExt.h"
#endif

template <class LStretchingInteractionType>
void LinkerStretching<LStretchingInteractionType>::assignforcemags() {
#ifdef CUDAACCL
    double stretchforce[Linker::getLinkers().size()];
    CUDAcommon::handleerror(cudaMemcpy(stretchforce, gpu_Lstretchforce,
                                       Linker::getLinkers().size() * sizeof(double),
                                       cudaMemcpyDeviceToHost));
    int id = 0;
    for(auto l:Linker::getLinkers())
    {l->getMLinker()->stretchForce = stretchforce[id];id++;}
#endif
}

template <class LStretchingInteractionType>
void LinkerStretching<LStretchingInteractionType>::vectorize() {

    beadSet = new int[n * Linker::getLinkers().size()];
    kstr = new double[Linker::getLinkers().size()];
    eql = new double[Linker::getLinkers().size()];
    pos1 = new double[Linker::getLinkers().size()];
    pos2 = new double[Linker::getLinkers().size()];
    stretchforce = new double[Linker::getLinkers().size()];

    int i = 0;

    for (auto l: Linker::getLinkers()) {
        l->_dbIndex = i;
        beadSet[n * i] = l->getFirstCylinder()->getFirstBead()->_dbIndex;
        beadSet[n * i + 1] = l->getFirstCylinder()->getSecondBead()->_dbIndex;
        beadSet[n * i + 2] = l->getSecondCylinder()->getFirstBead()->_dbIndex;
        beadSet[n * i + 3] = l->getSecondCylinder()->getSecondBead()->_dbIndex;

        kstr[i] = l->getMLinker()->getStretchingConstant();
        eql[i] = l->getMLinker()->getEqLength();
        pos1[i] = l->getFirstPosition();
        pos2[i] = l->getSecondPosition();
        stretchforce[i] = 0.0;
        i++;
    }

    //CUDA
#ifdef CUDAACCL
#ifdef CUDATIMETRACK
    chrono::high_resolution_clock::time_point tbegin, tend;
    tbegin = chrono::high_resolution_clock::now();
#endif
    //CUDA stream create
    if(stream == NULL || !(CUDAcommon::getCUDAvars().conservestreams))
        CUDAcommon::handleerror(cudaStreamCreate(&stream));
//    F_i = new double[3 * Bead::getBeads().size()];
//    cudaEvent_t start, stop;
//    CUDAcommon::handleerror(cudaEventCreate( &start));
//    CUDAcommon::handleerror(cudaEventCreate( &stop));
//    CUDAcommon::handleerror(cudaEventRecord( start, 0));

    int numInteractions =Linker::getLinkers().size();
    _FFType.optimalblocksnthreads(numInteractions, stream);
    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_beadSet, n * numInteractions * sizeof(int)),"cuda data transfer",
                                       "LinkerStretching.cu");
    CUDAcommon::handleerror(cudaMemcpyAsync(gpu_beadSet, beadSet, n * numInteractions *
                                                sizeof(int),
                                       cudaMemcpyHostToDevice, stream),"cuda data transfer", "LinkerStretching.cu");

    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_kstr, numInteractions * sizeof(double)),"cuda data transfer",
                                       "LinkerStretching.cu");
    CUDAcommon::handleerror(cudaMemcpyAsync(gpu_kstr, kstr, numInteractions * sizeof
                            (double), cudaMemcpyHostToDevice, stream),
                            "cuda data transfer", "LinkerStretching.cu");

    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_eql, numInteractions * sizeof(double)),"cuda data transfer",
                            "LinkerStretching.cu");
    CUDAcommon::handleerror(cudaMemcpyAsync(gpu_eql, eql, numInteractions * sizeof(double),
                            cudaMemcpyHostToDevice, stream), "cuda data transfer", "LinkerStretching.cu");
    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_pos1, numInteractions * sizeof(double)),"cuda data transfer",
                            "LinkerStretching.cu");
    CUDAcommon::handleerror(cudaMemcpyAsync(gpu_pos1, pos1, numInteractions * sizeof
                            (double), cudaMemcpyHostToDevice, stream),
                            "cuda data transfer", "LinkerStretching.cu");
    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_pos2, numInteractions * sizeof(double)),"cuda data transfer",
                            "LinkerStretching.cu");
    CUDAcommon::handleerror(cudaMemcpyAsync(gpu_pos2, pos2, numInteractions * sizeof
                            (double), cudaMemcpyHostToDevice, stream),
                            "cuda data transfer", "LinkerStretching.cu");
    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_Lstretchforce, numInteractions *
                                                                     sizeof(double)),"cuda data transfer",
                            "LinkerStretching.cu");
    vector<int> params;
    params.push_back(int(n));
    params.push_back(numInteractions);
    params.push_back(CUDAcommon::cudavars.offset_E);
    //set offset
    CUDAcommon::cudavars.offset_E += numInteractions;
//    std::cout<<"offset "<<getName()<<" "<<CUDAcommon::cudavars.offset_E<<endl;
    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_params, 3 * sizeof(int)),"cuda data"
                                    " transfer",
                            "LinkerStretching.cu");
    CUDAcommon::handleerror(cudaMemcpyAsync(gpu_params, params.data(), 3 * sizeof(int),
                                       cudaMemcpyHostToDevice, stream),
                            "cuda data transfer", "LinkerStretching.cu");
#ifdef CUDATIMETRACK
//    CUDAcommon::handleerror(cudaDeviceSynchronize(),"LinkerStretching.cu",
//                            "vectorizeFF");
    tend= chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed_run(tend - tbegin);
    CUDAcommon::cudatime.TvecvectorizeFF.push_back(elapsed_run.count());
    CUDAcommon::cudatime.TvectorizeFF += elapsed_run.count();
#endif
#endif
}

template<class LStretchingInteractionType>
void LinkerStretching<LStretchingInteractionType>::deallocate() {
    for(auto l:Linker::getLinkers()){
        //Using += to ensure that the stretching forces are additive.
        l->getMLinker()->stretchForce += stretchforce[l->_dbIndex];
    }
    delete [] stretchforce;
    delete [] beadSet;
    delete [] kstr;
    delete [] eql;
    delete [] pos1;
    delete [] pos2;
#ifdef CUDAACCL
        if(!(CUDAcommon::getCUDAvars().conservestreams))
            CUDAcommon::handleerror(cudaStreamDestroy(stream));
    _FFType.deallocate();
    CUDAcommon::handleerror(cudaFree(gpu_beadSet),"cudaFree", "LinkerStretching.cu");
    CUDAcommon::handleerror(cudaFree(gpu_kstr),"cudaFree", "LinkerStretching.cu");
    CUDAcommon::handleerror(cudaFree(gpu_pos1),"cudaFree", "LinkerStretching.cu");
    CUDAcommon::handleerror(cudaFree(gpu_pos2),"cudaFree", "LinkerStretching.cu");
    CUDAcommon::handleerror(cudaFree(gpu_eql),"cudaFree", "LinkerStretching.cu");
    CUDAcommon::handleerror(cudaFree(gpu_params),"cudaFree", "LinkerStretching.cu");
    CUDAcommon::handleerror(cudaFree(gpu_Lstretchforce),"cudaFree", "LinkerStretching.cu");
#endif
}


template <class LStretchingInteractionType>
double LinkerStretching<LStretchingInteractionType>::computeEnergy(double* coord, double *f, double d){

    double U_i[1], U_ii;
    double* gU_i;
    U_ii = 0.0;
#ifdef CUDATIMETRACK
    chrono::high_resolution_clock::time_point tbegin, tend;
#endif
#ifdef CUDAACCL
//    std::cout<<"Linker size "<<Linker::getLinkers().size()<<endl;
#ifdef CUDATIMETRACK
    tbegin = chrono::high_resolution_clock::now();
#endif

    //has to be changed to accomodate aux force
    double * gpu_coord=CUDAcommon::getCUDAvars().gpu_coord;
    double * gpu_force=CUDAcommon::getCUDAvars().gpu_force;
    double * gpu_d = CUDAcommon::getCUDAvars().gpu_lambda;


//    if(d == 0.0){
//        gU_i=_FFType.energy(gpu_coord, gpu_force, gpu_beadSet, gpu_kstr, gpu_eql, gpu_pos1, gpu_pos2, gpu_params);
//    }
//    else{
        gU_i=_FFType.energy(gpu_coord, gpu_force, gpu_beadSet, gpu_kstr, gpu_eql, gpu_pos1, gpu_pos2, gpu_d,
                            gpu_params);
//    }
#ifdef CUDATIMETRACK
//    CUDAcommon::handleerror(cudaDeviceSynchronize(),"CylinderExclVolume.cu",
//                            "computeEnergy");
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
        U_ii = _FFType.energy(coord, f, beadSet, kstr, eql, pos1, pos2);
    else
        U_ii = _FFType.energy(coord, f, beadSet, kstr, eql, pos1, pos2, d);

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

template <class LStretchingInteractionType>
void LinkerStretching<LStretchingInteractionType>::computeForces(double *coord, double *f) {
#ifdef CUDATIMETRACK
    chrono::high_resolution_clock::time_point tbegin, tend;
    tbegin = chrono::high_resolution_clock::now();
#endif
#ifdef CUDAACCL
    //has to be changed to accomodate aux force
    double * gpu_coord=CUDAcommon::getCUDAvars().gpu_coord;
    double * gpu_force;
    if(cross_checkclass::Aux){
        gpu_force=CUDAcommon::getCUDAvars().gpu_forceAux;
        _FFType.forces(gpu_coord, gpu_force, gpu_beadSet, gpu_kstr, gpu_eql, gpu_pos1,
                       gpu_pos2, gpu_params, gpu_Lstretchforce);
    }
    else {
        gpu_force = CUDAcommon::getCUDAvars().gpu_force;
        _FFType.forces(gpu_coord, gpu_force, gpu_beadSet, gpu_kstr, gpu_eql, gpu_pos1,
                       gpu_pos2, gpu_params , gpu_Lstretchforce);
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
    _FFType.forces(coord, f, beadSet, kstr, eql, pos1, pos2, stretchforce);
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
#ifdef CUDATIMETRACK
    tend= chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed_runs(tend - tbegin);
    CUDAcommon::serltime.TveccomputeF.push_back(elapsed_runs.count());
    CUDAcommon::serltime.TcomputeF += elapsed_runs.count();
#endif
}


///Temlate specializations
template double LinkerStretching<LinkerStretchingHarmonic>::computeEnergy(double *coord, double *f, double d);
template void LinkerStretching<LinkerStretchingHarmonic>::computeForces(double *coord, double *f);
template void LinkerStretching<LinkerStretchingHarmonic>::vectorize();
template void LinkerStretching<LinkerStretchingHarmonic>::deallocate();
template void LinkerStretching<LinkerStretchingHarmonic>::assignforcemags();
