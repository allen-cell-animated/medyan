
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

#include "MotorGhostStretching.h"

#include "MotorGhostStretchingHarmonic.h"

#include "MotorGhost.h"
#include "Cylinder.h"
#include "Bead.h"
#include "cross_check.h"
#include "CGMethod.h"
#ifdef CUDAACCL
#include "nvToolsExt.h"
#endif

template <class MStretchingInteractionType>
void MotorGhostStretching<MStretchingInteractionType>::assignforcemags() {
    // for(auto m: MotorGhost::getMotorGhosts()){
    //     //Using += to ensure that the stretching forces are additive.
    //     m->getMMotorGhost()->stretchForce = stretchforce[m->getIndex()];
    // }
#ifdef CUDAACCL
    floatingpoint stretchforce[MotorGhost::getMotorGhosts().size()];
    CUDAcommon::handleerror(cudaMemcpy(stretchforce, gpu_Mstretchforce,
                                       MotorGhost::getMotorGhosts().size() * sizeof(floatingpoint),
                                       cudaMemcpyDeviceToHost));
    int id = 0;
    for(auto m:MotorGhost::getMotorGhosts())
    {m->getMMotorGhost()->stretchForce = stretchforce[id];id++;}
#endif
}

template <class MStretchingInteractionType>
void MotorGhostStretching<MStretchingInteractionType>::vectorize(const FFCoordinateStartingIndex& si) {

    CUDAcommon::tmin.numinteractions[3] += MotorGhost::getMotorGhosts().size();
    beadSet = new int[n * MotorGhost::getMotorGhosts().size()];
    kstr = new floatingpoint[MotorGhost::getMotorGhosts().size()];
    eql = new floatingpoint[MotorGhost::getMotorGhosts().size()];
    pos1 = new floatingpoint[MotorGhost::getMotorGhosts().size()];
    pos2 = new floatingpoint[MotorGhost::getMotorGhosts().size()];
    stretchforce = new floatingpoint[MotorGhost::getMotorGhosts().size()];

    int i = 0;
    #ifdef ADDITIONALINFO
    MotorGhostInteractions::kstrvec.clear();
    MotorGhostInteractions::eqlvec.clear();
    #endif

    for (auto m: MotorGhost::getMotorGhosts()) {
        /* Haoran 03/18/2019 m->getIndex() = i; */
        beadSet[n * i] = m->getFirstCylinder()->getFirstBead()->getIndex() * 3 + si.bead;
        beadSet[n * i + 1] = m->getFirstCylinder()->getSecondBead()->getIndex() * 3 + si.bead;
        beadSet[n * i + 2] = m->getSecondCylinder()->getFirstBead()->getIndex() * 3 + si.bead;
        beadSet[n * i + 3] = m->getSecondCylinder()->getSecondBead()->getIndex() * 3 + si.bead;

        kstr[i] = m->getMMotorGhost()->getStretchingConstant();
        eql[i] = m->getMMotorGhost()->getEqLength();
        pos1[i] = m->getFirstCylinder()->adjustedrelativeposition(m->getFirstPosition());
        pos2[i] = m->getSecondCylinder()->adjustedrelativeposition(m->getSecondPosition());
        stretchforce[i] = 0.0;
        #ifdef ADDITIONALINFO
        MotorGhostInteractions::kstrvec.push_back(kstr[i]);
        MotorGhostInteractions::eqlvec.push_back(eql[i]);
		#endif

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
//    F_i = new floatingpoint[3 * Bead::getBeads().size()];
//    cudaEvent_t start, stop;
//    CUDAcommon::handleerror(cudaEventCreate( &start));
//    CUDAcommon::handleerror(cudaEventCreate( &stop));
//    CUDAcommon::handleerror(cudaEventRecord( start, 0));

    int numInteractions = MotorGhost::getMotorGhosts().size();
    _FFType.optimalblocksnthreads(numInteractions, stream);
//    blocksnthreads.clear();
//    blocksnthreads.push_back(numInteractions/THREADSPERBLOCK + 1);
//
//    if(blocksnthreads[0]==1) blocksnthreads.push_back( numInteractions);
////    if(blocksnthreads[0]==1) blocksnthreads.push_back( 32*(int(numInteractions/32 +1)) );
//    else blocksnthreads.push_back(THREADSPERBLOCK);
    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_beadSet, n * numInteractions * sizeof(int)));
    CUDAcommon::handleerror(cudaMemcpyAsync(gpu_beadSet, beadSet, n * numInteractions *
                                                sizeof(int),
                                       cudaMemcpyHostToDevice, stream));

    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_kstr, numInteractions * sizeof(floatingpoint)));
    CUDAcommon::handleerror(cudaMemcpyAsync(gpu_kstr, kstr, numInteractions * sizeof
                            (floatingpoint), cudaMemcpyHostToDevice, stream));

    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_eql, numInteractions * sizeof(floatingpoint)));
    CUDAcommon::handleerror(cudaMemcpyAsync(gpu_eql, eql, numInteractions * sizeof(floatingpoint),
                                        cudaMemcpyHostToDevice, stream));

    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_pos1, numInteractions * sizeof(floatingpoint)));
    CUDAcommon::handleerror(cudaMemcpyAsync(gpu_pos1, pos1, numInteractions * sizeof
                            (floatingpoint), cudaMemcpyHostToDevice, stream));

//    floatingpoint checkpos1[numInteractions];
//    cudaMemcpy(checkpos1, gpu_pos1, numInteractions * sizeof(floatingpoint), cudaMemcpyDeviceToHost);
//    for(auto i=0;i<numInteractions;i++) std::cout<<pos1[i]<<" "<<checkpos1[i]<<endl;

    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_pos2, numInteractions * sizeof(floatingpoint)));
    CUDAcommon::handleerror(cudaMemcpyAsync(gpu_pos2, pos2, numInteractions * sizeof
                           (floatingpoint), cudaMemcpyHostToDevice, stream));
    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_Mstretchforce, numInteractions *
                                                                     sizeof(floatingpoint)),"cuda data transfer",
                            "MotorGhostStretching.cu");

    vector<int> params;
    params.push_back(int(n));
    params.push_back(numInteractions);
    params.push_back(CUDAcommon::cudavars.offset_E);
    //set offset
    CUDAcommon::cudavars.offset_E += numInteractions;
//    std::cout<<"offset "<<getName()<<" "<<CUDAcommon::cudavars.offset_E<<endl;

    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_params, 3 * sizeof(int)));
    CUDAcommon::handleerror(cudaMemcpyAsync(gpu_params, params.data(), 3 * sizeof(int),
                                       cudaMemcpyHostToDevice, stream));
//    CUDAcommon::cudavars.motorparams = gpu_params;
#ifdef CUDATIMETRACK
//    CUDAcommon::handleerror(cudaDeviceSynchronize(),"MotorGhostStretching.cu",
//                            "vectorizeFF");
    tend= chrono::high_resolution_clock::now();
    chrono::duration<floatingpoint> elapsed_run(tend - tbegin);
    CUDAcommon::cudatime.TvecvectorizeFF.push_back(elapsed_run.count());
    CUDAcommon::cudatime.TvectorizeFF += elapsed_run.count();
#endif
#endif

    //
}

template<class MStretchingInteractionType>
void MotorGhostStretching<MStretchingInteractionType>::deallocate() {
    for(auto m: MotorGhost::getMotorGhosts()){
        //Using += to ensure that the stretching forces are additive.
        m->getMMotorGhost()->stretchForce += stretchforce[m->getIndex()];
//        std::cout<<m->getMMotorGhost()->stretchForce<<endl;
//        cout<<stretchforce[m->_dbIndex]<<" ";
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
    CUDAcommon::handleerror(cudaFree(gpu_beadSet));
    CUDAcommon::handleerror(cudaFree(gpu_kstr));
    CUDAcommon::handleerror(cudaFree(gpu_pos1));
    CUDAcommon::handleerror(cudaFree(gpu_pos2));
    CUDAcommon::handleerror(cudaFree(gpu_eql));
    CUDAcommon::handleerror(cudaFree(gpu_params));
    CUDAcommon::handleerror(cudaFree(gpu_Mstretchforce));
#endif
}


template <class MStretchingInteractionType>
floatingpoint MotorGhostStretching<MStretchingInteractionType>::computeEnergy(floatingpoint* coord){
    floatingpoint U_ii=0.0;
#ifdef CUDATIMETRACK
    chrono::high_resolution_clock::time_point tbegin, tend;
#endif
#ifdef CUDAACCL
//    std::cout<<"Motor size "<<MotorGhost::getMotorGhosts().size()<<endl;
#ifdef CUDATIMETRACK
    tbegin = chrono::high_resolution_clock::now();
#endif

    //has to be changed to accomodate aux force
    floatingpoint * gpu_coord=CUDAcommon::getCUDAvars().gpu_coord;
    floatingpoint * gpu_force=CUDAcommon::getCUDAvars().gpu_force;
    floatingpoint * gpu_d = CUDAcommon::getCUDAvars().gpu_lambda;


//    if(d == 0.0){
//        gU_i=_FFType.energy(gpu_coord, gpu_force, gpu_beadSet, gpu_kstr, gpu_eql, gpu_pos1, gpu_pos2, gpu_params);
//
//    }
//    else{
        gU_i=_FFType.energy(gpu_coord, gpu_force, gpu_beadSet, gpu_kstr, gpu_eql, gpu_pos1, gpu_pos2, gpu_d,
                            gpu_params);
//    }


#ifdef CUDATIMETRACK
//    CUDAcommon::handleerror(cudaDeviceSynchronize(),"MotorGhostStretching.cu", "computeEnergy");
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

    U_ii = _FFType.energy(coord, beadSet, kstr, eql, pos1, pos2);

#ifdef CUDATIMETRACK
    floatingpoint U_i[1], U_ii=0.0;
    floatingpoint* gU_i;
    tend= chrono::high_resolution_clock::now();
    chrono::duration<floatingpoint> elapsed_runs(tend - tbegin);
    CUDAcommon::serltime.TveccomputeE.push_back(elapsed_runs.count());
    CUDAcommon::serltime.TcomputeE += elapsed_runs.count();
    CUDAcommon::serltime.TcomputeEiter += elapsed_runs.count();
#endif
#endif

    return U_ii;
}

template <class MStretchingInteractionType>
void MotorGhostStretching<MStretchingInteractionType>::computeForces(floatingpoint *coord, floatingpoint *f) {
#ifdef CUDATIMETRACK
    chrono::high_resolution_clock::time_point tbegin, tend;
    tbegin = chrono::high_resolution_clock::now();
#endif
#ifdef CUDAACCL
    //has to be changed to accomodate aux force
    floatingpoint * gpu_coord=CUDAcommon::getCUDAvars().gpu_coord;

    floatingpoint * gpu_force;
    if(cross_checkclass::Aux){
        gpu_force=CUDAcommon::getCUDAvars().gpu_forceAux;
        _FFType.forces(gpu_coord, gpu_force, gpu_beadSet, gpu_kstr, gpu_eql, gpu_pos1,
                       gpu_pos2, gpu_params, gpu_Mstretchforce);
    }
    else {
        gpu_force = CUDAcommon::getCUDAvars().gpu_force;
        _FFType.forces(gpu_coord, gpu_force, gpu_beadSet, gpu_kstr, gpu_eql, gpu_pos1,
                       gpu_pos2, gpu_params, gpu_Mstretchforce);
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
    _FFType.forces(coord, f, beadSet, kstr, eql, pos1, pos2, stretchforce);
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


///Temlate specializations
template floatingpoint MotorGhostStretching<MotorGhostStretchingHarmonic>::computeEnergy(floatingpoint *coord);
template void MotorGhostStretching<MotorGhostStretchingHarmonic>::computeForces(floatingpoint *coord, floatingpoint *f);
template void MotorGhostStretching<MotorGhostStretchingHarmonic>::vectorize(const FFCoordinateStartingIndex&);
template void MotorGhostStretching<MotorGhostStretchingHarmonic>::deallocate();
template void MotorGhostStretching<MotorGhostStretchingHarmonic>::assignforcemags();


