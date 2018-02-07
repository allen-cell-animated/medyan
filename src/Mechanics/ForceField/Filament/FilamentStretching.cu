
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

#include "FilamentStretching.h"

#include "FilamentStretchingHarmonic.h"
#include "Bead.h"
#include "cross_check.h"
#include "nvToolsExt.h"

template <class FStretchingInteractionType>
void FilamentStretching<FStretchingInteractionType>::vectorize() {

    beadSet = new int[n * Cylinder::getCylinders().size()];
    kstr = new double[Cylinder::getCylinders().size()];
    eql = new double[Cylinder::getCylinders().size()];

    int i = 0;

    for (auto c: Cylinder::getCylinders()) {
        beadSet[n * i] = c->getFirstBead()->_dbIndex;
        beadSet[n * i + 1] = c->getSecondBead()->_dbIndex;

        kstr[i] = c->getMCylinder()->getStretchingConst();
        eql[i] = c->getMCylinder()->getEqLength();

        i++;
    }
    //CUDA
#ifdef CUDAACCL
    nvtxRangePushA("CVFF");
    F_i = new double[3 * Bead::getBeads().size()];
    int numInteractions = Cylinder::getCylinders().size();
    _FFType.optimalblocksnthreads(numInteractions);

    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_beadSet, n * numInteractions * sizeof(int)),"cuda data "
            "transfer", " FilamentStretching.cu");
    CUDAcommon::handleerror(cudaMemcpy(gpu_beadSet, beadSet, n * numInteractions * sizeof(int),
                                       cudaMemcpyHostToDevice),"cuda data transfer", " FilamentStretching.cu");

    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_kstr, numInteractions * sizeof(double)),"cuda data transfer",
                            " FilamentStretching.cu");
    CUDAcommon::handleerror(cudaMemcpy(gpu_kstr, kstr, numInteractions * sizeof(double), cudaMemcpyHostToDevice),
                            "cuda data transfer", " FilamentStretching.cu");

    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_eql, numInteractions * sizeof(double)),"cuda data transfer",
                            " FilamentStretching.cu");
    CUDAcommon::handleerror(cudaMemcpy(gpu_eql, eql, numInteractions * sizeof(double), cudaMemcpyHostToDevice),
                            "cuda data transfer", " FilamentStretching.cu");

    vector<int> params;
    params.push_back(int(n));
    params.push_back(numInteractions);
    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_params, 2 * sizeof(int)),"cuda data transfer",
                            " FilamentStretching.cu");
    CUDAcommon::handleerror(cudaMemcpy(gpu_params, params.data(), 2 * sizeof(int), cudaMemcpyHostToDevice),
            "cuda data transfer", " FilamentStretching.cu");
    nvtxRangePop();
#endif

    //
}

template<class FStretchingInteractionType>
void FilamentStretching<FStretchingInteractionType>::deallocate() {

    delete beadSet;
    delete kstr;
    delete eql;
#ifdef CUDAACCL
    _FFType.deallocate();
    CUDAcommon::handleerror(cudaFree(gpu_beadSet));
    CUDAcommon::handleerror(cudaFree(gpu_kstr));
    CUDAcommon::handleerror(cudaFree(gpu_eql));
    CUDAcommon::handleerror(cudaFree(gpu_params));
#endif
}


template <class FStretchingInteractionType>
double FilamentStretching<FStretchingInteractionType>::computeEnergy(double* coord, double *f, double d){

    double U_i[1], U_ii;
    double* gU_i;
    U_ii = NULL;
#ifdef CUDAACCL
    //has to be changed to accomodate aux force
    double * gpu_coord=CUDAcommon::getCUDAvars().gpu_coord;
    double * gpu_force=CUDAcommon::getCUDAvars().gpu_force;
    double * gpu_d = CUDAcommon::getCUDAvars().gpu_lambda;
    nvtxRangePushA("CCEFS");

    if(d == 0.0){
        gU_i=_FFType.energy(gpu_coord, gpu_force, gpu_beadSet, gpu_kstr, gpu_eql, gpu_params);

    }
    else{
        gU_i=_FFType.energy(gpu_coord, gpu_force, gpu_beadSet, gpu_kstr, gpu_eql, gpu_d,
                            gpu_params);
    }
    nvtxRangePop();
#endif
    nvtxRangePushA("SCEFS");

    if (d == 0.0)
        U_ii = _FFType.energy(coord, f, beadSet, kstr, eql);
    else
        U_ii = _FFType.energy(coord, f, beadSet, kstr, eql, d);
    nvtxRangePop();
    if(gU_i!=NULL) {

        CUDAcommon::handleerror(cudaMemcpy(U_i, gU_i, sizeof(double), cudaMemcpyDeviceToHost),
                                "computeEnergy", "FilamentStretching.cu");
    }
    else
        U_i[0] = 0.0;
    if(fabs(U_ii)>1000000.0) {
        if (fabs((U_ii - U_i[0]) / U_ii) > 0.0001){
            std::cout<<endl;
            std::cout << "CUDA FSE " << U_i[0] << endl;
            std::cout << "Vectorized FSE " << U_ii << endl;
            std::cout << "Precision match error" << fabs(U_ii - U_i[0]) << endl;
        }
    }
    else {
        if (fabs(U_ii - U_i[0]) > 1.0 / 100000000.0){
            std::cout << "Precision match " << fabs(U_ii - U_i[0]) << endl;
            std::cout << "CUDA FSE " << U_i[0] << endl;
            std::cout << "Vectorized FSE " << U_ii << endl;
//        exit(EXIT_FAILURE);
        }
    }

#ifdef CROSSCHECK
    double U2 = 0;
    double U_ii;
    for (auto f: Filament::getFilaments()) {

        U_ii = 0;

        if (d == 0.0){
            for(auto it : f->getCylinderVector()){

                Bead* b1 = it->getFirstBead();
                Bead* b2 = it->getSecondBead();
                double kStretch = it->getMCylinder()->getStretchingConst();
                double eqLength = it->getMCylinder()->getEqLength();
                U_ii += _FFType.energy(b1, b2, kStretch, eqLength);
            }
        }
        else {
            for(auto it : f->getCylinderVector()){
                Bead* b1 = it->getFirstBead();
                Bead* b2 = it->getSecondBead();

//                std::cout<<b1->coordinate[0]<<" "<<b1->coordinate[1]<<" "<<b1->coordinate[2]<<" "<<b2->coordinate[0]<<" "<<b2->coordinate[1]<<" "<<b2->coordinate[2]<<" ";
                double kStretch =it->getMCylinder()->getStretchingConst();
                double eqLength = it->getMCylinder()->getEqLength();

                U_ii += _FFType.energy(b1, b2, kStretch, eqLength, d);
            }
        }

        if(fabs(U_i) == numeric_limits<double>::infinity()
           || U_ii != U_ii || U_ii < -1.0) {

            //set culprit and return
            _filamentCulprit = f;

            U2=-1;
            break;
        }
        else
            U2 += U_ii;
    }
    if(abs(U_i-U2)<=U2/100000000000)
        std::cout<<"E S YES "<<endl;
    else
    {   std::cout<<U_i<<" "<<U2<<endl;
        exit(EXIT_FAILURE);
    }
#endif
    whoisCulprit();
    return U_ii;
}

template <class FStretchingInteractionType>
void FilamentStretching<FStretchingInteractionType>::computeForces(double *coord, double *f) {
#ifdef CUDAACCL
    //has to be changed to accomodate aux force
    double * gpu_coord=CUDAcommon::getCUDAvars().gpu_coord;

    double * gpu_force;


    if(cross_checkclass::Aux){
        nvtxRangePushA("CCFFS");

        gpu_force=CUDAcommon::getCUDAvars().gpu_forceAux;
        _FFType.forces(gpu_coord, gpu_force, gpu_beadSet, gpu_kstr, gpu_eql, gpu_params);
        nvtxRangePop();
    }
    else {
        nvtxRangePushA("CCFFS");

        gpu_force = CUDAcommon::getCUDAvars().gpu_force;
        _FFType.forces(gpu_coord, gpu_force, gpu_beadSet, gpu_kstr, gpu_eql, gpu_params);
        nvtxRangePop();
    }


    //TODO remove this later need not copy forces back to CPU.
    CUDAcommon::handleerror(cudaMemcpy(F_i, gpu_force, 3 * Bead::getBeads().size() *sizeof(double),
                                       cudaMemcpyDeviceToHost));
#endif
    nvtxRangePushA("SCFFS");
    _FFType.forces(coord, f, beadSet, kstr, eql);
    nvtxRangePop();
#ifdef CUDAACCL

    bool state = false;
    for(auto iter=0;iter<Bead::getBeads().size();iter++) {
        if (fabs(F_i[3 * iter] - f[3 * iter]) <=1.0/100000000.0 && fabs(F_i[3 * iter + 1] - f[3 * iter + 1])
                                                                   <=1.0/100000000.0 && fabs(F_i[3 * iter + 2] - f[3 * iter + 2]) <=1.0/100000000.0)
        {state = true;}
        else {
            state = false;
            std::cout<<endl;
            std::cout<<"FS Forces"<<endl;
            std::cout<<"Precision match "<<fabs(F_i[3 * iter] - f[3 * iter])<<" "<<fabs(F_i[3 * iter + 1] - f[3 *
                                                                                                              iter + 1])<<" "<<fabs(F_i[3 * iter + 2] - f[3 * iter + 2])<<endl;
            std::cout << "CUDA       " << F_i[3 * iter] << " " << F_i[3 * iter + 1] << " " << F_i[3 * iter + 2] << endl;
            std::cout << "Vectorized " << f[3 * iter] << " " << f[3 * iter + 1] << " " << f[3 * iter + 2] << endl;

//        exit(EXIT_FAILURE);
        }
    }
//    if(state)
//        std::cout<<"F M YES"<<endl;
#endif
#ifdef CROSSCHECK
    for (auto f: Filament::getFilaments()) {

        for(auto it : f->getCylinderVector()){

            Bead* b1 = it->getFirstBead();
            Bead* b2 = it->getSecondBead();
            double kStretch =it->getMCylinder()->getStretchingConst();
            double eqLength = it->getMCylinder()->getEqLength();


            if(cross_checkclass::Aux)
                _FFType.forcesAux(b1, b2, kStretch, eqLength);
            else
                _FFType.forces(b1, b2, kStretch, eqLength);


        }
    }
//    for(auto bd:Bead::getBeads())
//        std::cout<<bd->force[0]<<" "<<bd->force[1]<<" "<<bd->force[2]<<" ";
//    std::cout<<endl;

     if(cross_checkclass::Aux)
     {auto state=cross_check::crosscheckAuxforces(f);    std::cout<<"F S YES "<<state<<endl;}
     else
     { auto state=cross_check::crosscheckforces(f);     std::cout<<"F S YES "<<state<<endl;}

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
template double FilamentStretching<FilamentStretchingHarmonic>::computeEnergy(double *coord, double *f, double d);
template void FilamentStretching<FilamentStretchingHarmonic>::computeForces(double *coord, double *f);
template void FilamentStretching<FilamentStretchingHarmonic>::vectorize();
template void FilamentStretching<FilamentStretchingHarmonic>::deallocate();
//template void FilamentStretching<FilamentStretchingHarmonic>::whoisCulprit();
