
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

#include "FilamentStretching.h"

#include "FilamentStretchingHarmonic.h"
#include "Bead.h"
#include "cross_check.h"
#include "CGMethod.h"
#ifdef CUDAACCL
#include "nvToolsExt.h"
#endif

template <class FStretchingInteractionType>
void FilamentStretching<FStretchingInteractionType>::vectorize(const FFCoordinateStartingIndex& si) {
    CUDAcommon::tmin.numinteractions[0] += Cylinder::getCylinders().size();
    beadSet = new int[n * Cylinder::getCylinders().size()];
    kstr = new floatingpoint[Cylinder::getCylinders().size()];
    eql = new floatingpoint[Cylinder::getCylinders().size()];

    int i = 0;

    for (auto c: Cylinder::getCylinders()) {
        beadSet[n * i] = c->getFirstBead()->getIndex() * 3 + si.bead;
        beadSet[n * i + 1] = c->getSecondBead()->getIndex() * 3 + si.bead;
        kstr[i] = c->getMCylinder()->getStretchingConst();
        eql[i] = c->getMCylinder()->getEqLength();
/*        std::cout<<"Filstretching with cindex "<<c->_dcIndex<<" and ID "
                ""<<c->getID()<<" with bindices "<<c->getFirstBead()
                         ->getIndex()<<" "<<c->getSecondBead()->getIndex()<<endl;*/
        i++;
    }
    //CUDA
#ifdef CUDAACCL
#ifdef CUDATIMETRACK
    //stream create
    if(stream == NULL || !(CUDAcommon::getCUDAvars().conservestreams))
        CUDAcommon::handleerror(cudaStreamCreate(&stream));

    chrono::high_resolution_clock::time_point tbegin, tend;
    tbegin = chrono::high_resolution_clock::now();
#endif
    int numInteractions = Cylinder::getCylinders().size();
    _FFType.optimalblocksnthreads(numInteractions, stream);

    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_beadSet, n * numInteractions * sizeof(int)),"cuda data "
            "transfer", " FilamentStretching.cu");
    CUDAcommon::handleerror(cudaMemcpyAsync(gpu_beadSet, beadSet, n * numInteractions *
                            sizeof(int), cudaMemcpyHostToDevice, stream),"cuda data transfer", " FilamentStretching.cu");

    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_kstr, numInteractions * sizeof(floatingpoint)),"cuda data transfer",
                            " FilamentStretching.cu");
    CUDAcommon::handleerror(cudaMemcpyAsync(gpu_kstr, kstr, numInteractions * sizeof
                            (floatingpoint), cudaMemcpyHostToDevice, stream), "cuda data transfer", " FilamentStretching.cu");

    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_eql, numInteractions * sizeof(floatingpoint)),"cuda data transfer",
                            " FilamentStretching.cu");
    CUDAcommon::handleerror(cudaMemcpyAsync(gpu_eql, eql, numInteractions * sizeof(floatingpoint),
                            cudaMemcpyHostToDevice, stream), "cuda data transfer", " FilamentStretching.cu");

    vector<int> params;
    params.push_back(int(n));
    params.push_back(numInteractions);
    params.push_back(CUDAcommon::cudavars.offset_E);
    //set offset
    CUDAcommon::cudavars.offset_E += numInteractions;
//    std::cout<<"offset "<<getName()<<" "<<CUDAcommon::cudavars.offset_E<<endl;
    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_params, 3 * sizeof(int)),"cuda data"
                                    " transfer", " FilamentStretching.cu");
    CUDAcommon::handleerror(cudaMemcpyAsync(gpu_params, params.data(), 3 * sizeof(int),
                                       cudaMemcpyHostToDevice, stream),
            "cuda data transfer", " FilamentStretching.cu");
#ifdef CUDATIMETRACK
//    CUDAcommon::handleerror(cudaDeviceSynchronize(),"FilamentStretching.cu",
//                            "vectorizeFF");
    tend= chrono::high_resolution_clock::now();
    chrono::duration<floatingpoint> elapsed_run(tend - tbegin);
    CUDAcommon::cudatime.TvecvectorizeFF.push_back(elapsed_run.count());
    CUDAcommon::cudatime.TvectorizeFF += elapsed_run.count();
#endif
#endif

    //
}

template<class FStretchingInteractionType>
void FilamentStretching<FStretchingInteractionType>::deallocate() {

    delete [] beadSet;
    delete [] kstr;
    delete [] eql;
#ifdef CUDAACCL
    if(!(CUDAcommon::getCUDAvars().conservestreams))
        CUDAcommon::handleerror(cudaStreamDestroy(stream));
    _FFType.deallocate();
    CUDAcommon::handleerror(cudaFree(gpu_beadSet));
    CUDAcommon::handleerror(cudaFree(gpu_kstr));
    CUDAcommon::handleerror(cudaFree(gpu_eql));
    CUDAcommon::handleerror(cudaFree(gpu_params));
#endif
}


template <class FStretchingInteractionType>
floatingpoint FilamentStretching<FStretchingInteractionType>::computeEnergy(floatingpoint* coord){

    floatingpoint U_i[1], U_ii;
    floatingpoint* gU_i;
    U_ii=0.0;
#ifdef CUDAACCL
#ifdef CUDATIMETRACK
    chrono::high_resolution_clock::time_point tbegin, tend;
    tbegin = chrono::high_resolution_clock::now();
#endif
    //has to be changed to accomodate aux force
    floatingpoint * gpu_coord=CUDAcommon::getCUDAvars().gpu_coord;
    floatingpoint * gpu_force=CUDAcommon::getCUDAvars().gpu_force;
    floatingpoint * gpu_d = CUDAcommon::getCUDAvars().gpu_lambda;
//    std::cout<<"Fil Stretching Forces"<<endl;


//    if(d == 0.0){
//        gU_i=_FFType.energy(gpu_coord, gpu_force, gpu_beadSet, gpu_kstr, gpu_eql, gpu_params);
//
//    }
//    else{
        gU_i=_FFType.energy(gpu_coord, gpu_force, gpu_beadSet, gpu_kstr, gpu_eql, gpu_d,
                            gpu_params);
//    }
#ifdef CUDATIMETRACK
//    CUDAcommon::handleerror(cudaDeviceSynchronize(),"CylinderExclVolume.cu",
//                            "computeEnergy");
    tend= chrono::high_resolution_clock::now();
    chrono::duration<floatingpoint> elapsed_run(tend - tbegin);
    CUDAcommon::cudatime.TveccomputeE.push_back(elapsed_run.count());
    CUDAcommon::cudatime.TcomputeE += elapsed_run.count();
    CUDAcommon::cudatime.TcomputeEiter += elapsed_run.count();
#endif
#endif
#ifdef SERIAL
#ifdef CUDATIMETRACK
    chrono::high_resolution_clock::time_point tbegin, tend;
    tbegin = chrono::high_resolution_clock::now();
#endif

    U_ii = _FFType.energy(coord, beadSet, kstr, eql);

#ifdef CUDATIMETRACK
    tend = chrono::high_resolution_clock::now();
    chrono::duration<floatingpoint> elapsed_runs(tend - tbegin);
    CUDAcommon::serltime.TveccomputeE.push_back(elapsed_runs.count());
    CUDAcommon::serltime.TcomputeE += elapsed_runs.count();
    CUDAcommon::serltime.TcomputeEiter += elapsed_runs.count();
#endif
#endif
//    whoisCulprit();
    return U_ii;
}

template <class FStretchingInteractionType>
void FilamentStretching<FStretchingInteractionType>::computeForces(floatingpoint *coord, floatingpoint *f) {
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
        _FFType.forces(gpu_coord, gpu_force, gpu_beadSet, gpu_kstr, gpu_eql, gpu_params);
    }
    else {
        gpu_force = CUDAcommon::getCUDAvars().gpu_force;
        _FFType.forces(gpu_coord, gpu_force, gpu_beadSet, gpu_kstr, gpu_eql, gpu_params);
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
    _FFType.forces(coord, f, beadSet, kstr, eql);
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
    std::cout<<"max "<<getName()<<" "<<maxF<<" nint "<<Cylinder::getCylinders().size()
             <<endl;
#endif
}

//template <class FStretchingInteractionType>
//void FilamentStretching<FStretchingInteractionType>::whoisCulprit() {
//    cudaDeviceSynchronize();
//    std::cout<<cudaGetLastError()<<endl;
//    cout<<cudaGetErrorString(cudaGetLastError())<<endl;
//    if( cudaGetLastError() == cudaErrorAssert){
//        cudaDeviceSynchronize();
//        auto culpritinteraction = CUDAcommon::getCUDAvars().culpritinteraction;
//        std::cout<<strcmp(culpritinteraction,"Filament Stretching")<<endl;
//        if(strcmp(culpritinteraction, "Filament Stretching")==0){
//            CUDAcommon::printculprit("FilamentStretching","FilamentStretchingHarmonic");
//            Filament* f;
//            f = (Filament*)(Cylinder::getCylinders()[CUDAcommon::getCUDAvars().culpritID[0]]->getParent());
//            cout<<"Printing culprit Filament information."<<endl;
//            f->printSelf();
//            exit(EXIT_FAILURE);
//        }
//    }
//}
///Temlate specializations
template floatingpoint FilamentStretching<FilamentStretchingHarmonic>::computeEnergy(floatingpoint *coord);
template void FilamentStretching<FilamentStretchingHarmonic>::computeForces(floatingpoint *coord, floatingpoint *f);
template void FilamentStretching<FilamentStretchingHarmonic>::vectorize();
template void FilamentStretching<FilamentStretchingHarmonic>::deallocate();
//template void FilamentStretching<FilamentStretchingHarmonic>::whoisCulprit();
