
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

#include "MotorGhostStretching.h"

#include "MotorGhostStretchingHarmonic.h"

#include "MotorGhost.h"
#include "Cylinder.h"
#include "Bead.h"
#include "cross_check.h"

template <class MStretchingInteractionType>
void MotorGhostStretching<MStretchingInteractionType>::vectorize() {

    beadSet = new int[n * MotorGhost::getMotorGhosts().size()];
    kstr = new double[MotorGhost::getMotorGhosts().size()];
    eql = new double[MotorGhost::getMotorGhosts().size()];
    pos1 = new double[MotorGhost::getMotorGhosts().size()];
    pos2 = new double[MotorGhost::getMotorGhosts().size()];

    int i = 0;

    for (auto m: MotorGhost::getMotorGhosts()) {

        beadSet[n * i] = m->getFirstCylinder()->getFirstBead()->_dbIndex;
        beadSet[n * i + 1] = m->getFirstCylinder()->getSecondBead()->_dbIndex;
        beadSet[n * i + 2] = m->getSecondCylinder()->getFirstBead()->_dbIndex;
        beadSet[n * i + 3] = m->getSecondCylinder()->getSecondBead()->_dbIndex;

        kstr[i] = m->getMMotorGhost()->getStretchingConstant();
        eql[i] = m->getMMotorGhost()->getEqLength();
        pos1[i] = m->getFirstPosition();
        pos2[i] = m->getSecondPosition();

        i++;
    }


    //CUDA
#ifdef CUDAACCL
    int numInteractions = MotorGhost::getMotorGhosts().size();
    blocksnthreads.clear();
    blocksnthreads.push_back(numInteractions/THREADSPERBLOCK + 1);

    if(blocksnthreads[0]==1) blocksnthreads.push_back( numInteractions);
//    if(blocksnthreads[0]==1) blocksnthreads.push_back( 32*(int(numInteractions/32 +1)) );
    else blocksnthreads.push_back(THREADSPERBLOCK);
    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_beadSet, n * numInteractions * sizeof(int)));
    CUDAcommon::handleerror(cudaMemcpy(gpu_beadSet, beadSet, n * numInteractions * sizeof(int),
                                       cudaMemcpyHostToDevice));

    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_kstr, numInteractions * sizeof(double)));
    CUDAcommon::handleerror(cudaMemcpy(gpu_kstr, kstr, numInteractions * sizeof(double), cudaMemcpyHostToDevice));

    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_eql, numInteractions * sizeof(double)));
    CUDAcommon::handleerror(cudaMemcpy(gpu_eql, eql, numInteractions * sizeof(double), cudaMemcpyHostToDevice));

    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_pos1, numInteractions * sizeof(double)));
    CUDAcommon::handleerror(cudaMemcpy(gpu_pos1, pos1, numInteractions * sizeof(double), cudaMemcpyHostToDevice));

//    double checkpos1[numInteractions];
//    cudaMemcpy(checkpos1, gpu_pos1, numInteractions * sizeof(double), cudaMemcpyDeviceToHost);
//    for(auto i=0;i<numInteractions;i++) std::cout<<pos1[i]<<" "<<checkpos1[i]<<endl;

    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_pos2, numInteractions * sizeof(double)));
    CUDAcommon::handleerror(cudaMemcpy(gpu_pos2, pos2, numInteractions * sizeof(double), cudaMemcpyHostToDevice));

    vector<int> params;
    params.push_back(int(n));
    params.push_back(numInteractions);
    std::cout<<params[0]<<" "<<params[1]<<endl;

    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_params, 2 * sizeof(int)));
    CUDAcommon::handleerror(cudaMemcpy(gpu_params, params.data(), 2 * sizeof(int), cudaMemcpyHostToDevice));
    CUDAcommon::cudavars.motorparams = gpu_params;

#endif
    //
}

template<class MStretchingInteractionType>
void MotorGhostStretching<MStretchingInteractionType>::deallocate() {

    delete beadSet;
    delete kstr;
    delete eql;
    delete pos1;
    delete pos2;
#ifdef CUDAACCL
    CUDAcommon::handleerror(cudaFree(gpu_beadSet));
    CUDAcommon::handleerror(cudaFree(gpu_kstr));
    CUDAcommon::handleerror(cudaFree(gpu_pos1));
    CUDAcommon::handleerror(cudaFree(gpu_pos2));
    CUDAcommon::handleerror(cudaFree(gpu_eql));
    CUDAcommon::handleerror(cudaFree(gpu_params));
#endif
}


template <class MStretchingInteractionType>
double MotorGhostStretching<MStretchingInteractionType>::computeEnergy(double* coord, double *f, double d){
    double U_i, U_ii;
#ifdef CUDAACCL
    //has to be changed to accomodate aux force
    double * gpu_coord=CUDAcommon::getCUDAvars().gpu_coord;
    double * gpu_force=CUDAcommon::getCUDAvars().gpu_force;
    double * gpu_d = CUDAcommon::getCUDAvars().gpu_lambda;
    std::cout<<"cuda----"<<endl;
    if(d == 0.0){

        //CROSSCHECK
//        std::cout<<"Number of beads "<<Bead::getBeads().size()<<endl;
//        for(auto i=0;i<3 * Bead::getBeads().size();i++){
//            std::cout<<ccoord[i]<<" ";
//            if(i%3==2) std::cout<<endl;
//        }
//        std::cout<<"---------END-----"<<endl;
        //
        U_i=_FFType.energy(gpu_coord, gpu_force, gpu_beadSet, gpu_kstr, gpu_eql, gpu_pos1, gpu_pos2, gpu_params, blocksnthreads);

    }
    else{
//        double dd[1];
//        CUDAcommon::handleerror(cudaMemcpy(dd, gpu_d, sizeof(double), cudaMemcpyDeviceToHost));
//        std::cout<<"d = "<<dd[0]<<" "<<d<<endl;
//        int cparams[2];
//        CUDAcommon::handleerror(cudaMemcpy(cparams, gpu_params, 2*sizeof(int), cudaMemcpyDeviceToHost));
//        std::cout<<"cparams "<<cparams[0]<<" "<<cparams[1]<<endl;


        U_i=_FFType.energy(gpu_coord, gpu_force, gpu_beadSet, gpu_kstr, gpu_eql, gpu_pos1, gpu_pos2, gpu_d, gpu_params,
                           blocksnthreads);
    }
#endif
    std::cout<<"vector----"<<endl;
    if (d == 0.0)
        U_ii = _FFType.energy(coord, f, beadSet, kstr, eql, pos1, pos2);
    else
        U_ii = _FFType.energy(coord, f, beadSet, kstr, eql, pos1, pos2, d);

    if(fabs(U_ii)>1000000.0) {
        if (fabs((U_ii - U_i) / U_ii) < 0.0001)
            std::cout << "E V+M YES" << endl;
        else {
            std::cout << "CUDA MSE " << U_i << endl;
            std::cout << "Vectorized MSE " << U_ii << endl;
            std::cout << "Precision match %error" << fabs((U_ii - U_i) / U_ii) << endl;
        }
    }
    else {
        if (fabs(U_ii - U_i) < 1.0 / 100000000.0)
            std::cout << "E V+M YES" << endl;
        else {
            std::cout << "CUDA MSE " << U_i << endl;
            std::cout << "Vectorized MSE " << U_ii << endl;
            std::cout << "Precision match " << fabs(U_ii - U_i) << endl;
//        exit(EXIT_FAILURE);
        }
    }

#ifdef CROSSCHECK
    double U2 = 0;
    double U_ii;

    for (auto m: MotorGhost::getMotorGhosts()) {

        Bead* b1 = m->getFirstCylinder()->getFirstBead();
        Bead* b2 = m->getFirstCylinder()->getSecondBead();
        Bead* b3 = m->getSecondCylinder()->getFirstBead();
        Bead* b4 = m->getSecondCylinder()->getSecondBead();

        double kStretch = m->getMMotorGhost()->getStretchingConstant();
        double eqLength = m->getMMotorGhost()->getEqLength();

        double pos1 = m->getFirstPosition();
        double pos2 = m->getSecondPosition();

        if (d == 0.0)
            U_ii = _FFType.energy(b1, b2, b3, b4, pos1, pos2, kStretch, eqLength);
        else
            U_ii = _FFType.energy(b1, b2, b3, b4, pos1, pos2, kStretch, eqLength, d);

        if(fabs(U_i) == numeric_limits<double>::infinity()
           || U_ii != U_ii || U_ii < -1.0) {

            //set culprit and return
            _motorCulprit = m;

            U2=-1;
            break;
        }
        else
            U2 += U_ii;
    }
    if(abs(U_i-U2)<=U2/100000000000)
        std::cout<<"E M YES "<<endl;
    else
    {   std::cout<<U_i<<" "<<U2<<endl;
        exit(EXIT_FAILURE);
    }


#endif

    return U_i;
}

template <class MStretchingInteractionType>
void MotorGhostStretching<MStretchingInteractionType>::computeForces(double *coord, double *f) {

#ifdef CUDAACCL
    //has to be changed to accomodate aux force
    double * gpu_coord=CUDAcommon::getCUDAvars().gpu_coord;

    double * gpu_force;

//    if(cross_checkclass::Aux) {
//        gpu_force = CUDAcommon::getCUDAvars().gpu_forceAux;
//    }
//    else
//        gpu_force = CUDAcommon::getCUDAvars().gpu_force;
//    double F_c[3*Bead::getBeads().size()];
//    //TODO remove this later need not copy forces back to CPU.
//    CUDAcommon::handleerror(cudaMemcpy(F_c, gpu_force, 3 * Bead::getBeads().size() *sizeof(double),
//                                       cudaMemcpyDeviceToHost));
//    cout.precision(dbl::max_digits10);
//    for(int iter=0;iter<Bead::getBeads().size();iter++) {
//        std::cout << "C " << F_c[3 * iter] << " " << F_c[3 * iter + 1] << " " << F_c[3 * iter + 2] <<" ";
//        std::cout << "V "<<f[3 * iter] << " " << f[3 * iter + 1] << " " << f[3 * iter + 2] << endl;
//    }
//    std::cout<<"check ends"<<endl;

    if(cross_checkclass::Aux){
        gpu_force=CUDAcommon::getCUDAvars().gpu_forceAux;
        _FFType.forces(gpu_coord, gpu_force, gpu_beadSet, gpu_kstr, gpu_eql, gpu_pos1, gpu_pos2, gpu_params,
                       blocksnthreads);
    }
    else {
        gpu_force = CUDAcommon::getCUDAvars().gpu_force;
        _FFType.forces(gpu_coord, gpu_force, gpu_beadSet, gpu_kstr, gpu_eql, gpu_pos1, gpu_pos2, gpu_params,
                       blocksnthreads);
    }

    double F_i[3*Bead::getBeads().size()];
    //TODO remove this later need not copy forces back to CPU.
    CUDAcommon::handleerror(cudaMemcpy(F_i, gpu_force, 3 * Bead::getBeads().size() *sizeof(double),
                                       cudaMemcpyDeviceToHost));
#endif
    std::cout<<"vectorized "<<endl;
    _FFType.forces(coord, f, beadSet, kstr, eql, pos1, pos2);

#ifdef CUDAACCL
    cout.precision(dbl::max_digits10);
//    std::cout<<"M forces"<<endl;
//    for(int iter=0;iter<Bead::getBeads().size();iter++) {
//        std::cout << F_i[3 * iter] << " " << F_i[3 * iter + 1] << " " << F_i[3 * iter + 2] <<" ";
//        std::cout <<f[3 * iter] << " " << f[3 * iter + 1] << " " << f[3 * iter + 2] << endl;
//    }

    bool state = false;
    for(int iter=0;iter<Bead::getBeads().size();iter++) {
        if (fabs(F_i[3 * iter] - f[3 * iter]) <=1.0/100000000.0 && fabs(F_i[3 * iter + 1] - f[3 * iter + 1])
        <=1.0/100000000.0 && fabs(F_i[3 * iter + 2] - f[3 * iter + 2]) <=1.0/100000000.0)
        {state = true;}
        else {
//            std::cout<<endl;
            std::cout << "CUDA       " << F_i[3 * iter] << " " << F_i[3 * iter + 1] << " " << F_i[3 * iter + 2] << endl;
            std::cout << "Vectorized " << f[3 * iter] << " " << f[3 * iter + 1] << " " << f[3 * iter + 2] << endl;
            std::cout<<"Precision match "<<fabs(F_i[3 * iter] - f[3 * iter])<<" "<<fabs(F_i[3 * iter + 1] - f[3 *
            iter + 1])<<" "<<fabs(F_i[3 * iter + 2] - f[3 * iter + 2])<<endl;
//        exit(EXIT_FAILURE);
        }
    }
    if(state)
        std::cout<<endl;
        std::cout<<"F M YES"<<endl;
#endif

#ifdef CROSSCHECK
    for (auto m: MotorGhost::getMotorGhosts()) {

        Bead* b1 = m->getFirstCylinder()->getFirstBead();
        Bead* b2 = m->getFirstCylinder()->getSecondBead();
        Bead* b3 = m->getSecondCylinder()->getFirstBead();
        Bead* b4 = m->getSecondCylinder()->getSecondBead();
        double kStretch = m->getMMotorGhost()->getStretchingConstant();
        double eqLength = m->getMMotorGhost()->getEqLength();

        double pos1 = m->getFirstPosition();
        double pos2 = m->getSecondPosition();

        if(cross_checkclass::Aux)
        {double f0 = _FFType.forcesAux(b1, b2, b3, b4, pos1, pos2, kStretch, eqLength);
//            m->getMMotorGhost()->stretchForce = f0;
        }
        else
        {double f0 = _FFType.forces(b1, b2, b3, b4, pos1, pos2, kStretch, eqLength);
//            m->getMMotorGhost()->stretchForce = f0;
        }
    }
    if(cross_checkclass::Aux){
        auto state=cross_check::crosscheckAuxforces(f);
        std::cout<<"F S+B+L+M YES "<<state<<endl;}
    else{
        auto state=cross_check::crosscheckforces(f);
        std::cout<<"F S+B+L+M YES "<<state<<endl;}
#endif
}


///Temlate specializations
template double MotorGhostStretching<MotorGhostStretchingHarmonic>::computeEnergy(double *coord, double *f, double d);
template void MotorGhostStretching<MotorGhostStretchingHarmonic>::computeForces(double *coord, double *f);
template void MotorGhostStretching<MotorGhostStretchingHarmonic>::vectorize();
template void MotorGhostStretching<MotorGhostStretchingHarmonic>::deallocate();


