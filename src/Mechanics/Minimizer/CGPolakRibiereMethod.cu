

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

#include "CGPolakRibiereMethod.h"

#include "ForceFieldManager.h"
#include "Composite.h"
#include "Output.h"
#include "cross_check.h"
void PolakRibiere::minimize(ForceFieldManager &FFM, double GRADTOL,
                            double MAXDIST, double LAMBDAMAX, bool steplimit){

    //number of steps
    int N;
    if(steplimit) {
        int beadMaxStep = 5 * Bead::numBeads();
        N = (beadMaxStep > _MINNUMSTEPS ? beadMaxStep : _MINNUMSTEPS);
    }
    else {
        N = numeric_limits<int>::max();
    }

    startMinimization();
//    long index = 0;
//    long i = 0;
//    for(auto b: Bead::getBeads()) {
//
//        //flatten indices
//        index = 3 * i;
//        std::cout<<coord[index]<<" "<<coord[index+1]<<" "<<coord[index+2]<<" ";
//
//        i++;
//    }
//    std::cout<<endl;

    FFM.vectorizeAllForceFields();
#if defined(CROSSCHECK) || defined(CUDAENABLED) ||defined(CUDACROSSCHECK)
    cross_checkclass::Aux=false;
#endif
    FFM.computeForces(coord, force);
    FFM.copyForces(forceAux, force);

    double F_i[3*Bead::getBeads().size()];
    double max2;
    CUDAcommon::handleerror(cudaMemcpy(F_i, CUDAcommon::getCUDAvars().gpu_force, 3 * Bead::getBeads().size() *
                                                                                 sizeof(double), cudaMemcpyDeviceToHost));
    max2=0.0;
    std::cout<<"forces"<<endl;
    for(int iter=0;iter<Bead::getBeads().size();iter++) {
//        std::cout << "C " << F_i[3 * iter] << " " << F_i[3 * iter + 1] << " " << F_i[3 * iter + 2] <<" ";
//        std::cout << "V "<<force[3 * iter] << " " << force[3 * iter + 1] << " " << force[3 * iter + 2] << endl;
        if(abs(F_i[3 * iter]-force[3 * iter])>max2)
            max2 = abs(F_i[3 * iter]-force[3 * iter]);
        else if(abs(F_i[3 * iter +1]-force[3 * iter +1])>max2)
            max2 = abs(F_i[3 * iter +1]-force[3 * iter +1]);
        else if(abs(F_i[3 * iter +2]-force[3 * iter +2])>max2)
            max2 = abs(F_i[3 * iter +1]-force[3 * iter +1]);
    }

    std::cout<<"force alone start "<<max2<<endl;
    //TODO write the following function in CUDA
    double curGrad = CGMethod::allFDotF();

    int numIter = 0;

    while (/* Iteration criterion */  numIter < N &&
                                      /* Gradient tolerance  */  maxF() > GRADTOL) {
        std::cout<<"CONDITION "<< numIter<<" "<<maxF()<<" "<<GRADTOL<<endl;
        numIter++;
        double lambda, beta, newGrad, prevGrad;

#if defined(CROSSCHECK)
        auto state=cross_check::crosscheckforces(force);
#endif

        //find lambda by line search, move beads
        lambda = _safeMode ? safeBacktrackingLineSearch(FFM, MAXDIST, LAMBDAMAX)
                           : backtrackingLineSearch(FFM, MAXDIST, LAMBDAMAX);

#ifdef CUDAACCL

        // TODO change endif to else so that moveBeads affects gpu_coord.
#endif
        double F_i[3*Bead::getBeads().size()];
        double C_i[3*Bead::getBeads().size()];
        double max1,max2,max3;
        CUDAcommon::handleerror(cudaMemcpy(C_i, CUDAcommon::getCUDAvars().gpu_coord, 3 * Bead::getBeads().size() *
                                                                                     sizeof(double), cudaMemcpyDeviceToHost));
        max1=0.0;
        for(int iter=0;iter<Bead::getBeads().size();iter++) {
//            std::cout << "C " << C_i[3 * iter] << " " << C_i[3 * iter + 1] << " " << C_i[3 * iter + 2] <<" ";
//            std::cout << "V "<<coord[3 * iter] << " " << coord[3 * iter + 1] << " " <<coord[3 * iter + 2] << endl;
            if(abs(C_i[3 * iter]-coord[3 * iter])>max1)
                max1 = abs(C_i[3 * iter]-coord[3 * iter]);
            else if(abs(C_i[3 * iter +1]-coord[3 * iter +1])>max1)
                max1 = abs(C_i[3 * iter +1]-coord[3 * iter +1]);
            else if(abs(C_i[3 * iter +2]-coord[3 * iter +2])>max1)
                max1 = abs(C_i[3 * iter +1]-coord[3 * iter +1]);
        }

        moveBeads(lambda);

        CUDAcommon::handleerror(cudaMemcpy(C_i, CUDAcommon::getCUDAvars().gpu_coord, 3 * Bead::getBeads().size() *
                                                                                     sizeof(double), cudaMemcpyDeviceToHost));
        CUDAcommon::handleerror(cudaMemcpy(F_i, CUDAcommon::getCUDAvars().gpu_force, 3 * Bead::getBeads().size() *
                                                                                     sizeof(double), cudaMemcpyDeviceToHost));
        max2=0.0;max3=0.0;
        std::cout<<"forces"<<endl;
        for(int iter=0;iter<Bead::getBeads().size();iter++) {
//            std::cout << "C " << F_i[3 * iter] << " " << F_i[3 * iter + 1] << " " << F_i[3 * iter + 2] <<" ";
//            std::cout << "V "<<force[3 * iter] << " " << force[3 * iter + 1] << " " << force[3 * iter + 2] << endl;
            if(abs(F_i[3 * iter]-force[3 * iter])>max2)
                max2 = abs(F_i[3 * iter]-force[3 * iter]);
            else if(abs(F_i[3 * iter +1]-force[3 * iter +1])>max2)
                max2 = abs(F_i[3 * iter +1]-force[3 * iter +1]);
            else if(abs(F_i[3 * iter +2]-force[3 * iter +2])>max2)
                max2 = abs(F_i[3 * iter +1]-force[3 * iter +1]);
        }
        std::cout<<"coord"<<endl;
        for(int iter=0;iter<Bead::getBeads().size();iter++) {
//            std::cout << "C " << C_i[3 * iter] << " " << C_i[3 * iter + 1] << " " << C_i[3 * iter + 2] <<" ";
//            std::cout << "V "<<coord[3 * iter] << " " << coord[3 * iter + 1] << " " <<coord[3 * iter + 2] << endl;
            if(abs(C_i[3 * iter]-coord[3 * iter])>max3)
                max3 = abs(C_i[3 * iter]-coord[3 * iter]);
            else if(abs(C_i[3 * iter +1]-coord[3 * iter +1])>max3)
                max3 = abs(C_i[3 * iter +1]-coord[3 * iter +1]);
            else if(abs(C_i[3 * iter +2]-coord[3 * iter +2])>max3)
                max3 = abs(C_i[3 * iter +1]-coord[3 * iter +1]);
        }
        double gchecklambda[1];
        CUDAcommon::handleerror(cudaMemcpy(gchecklambda, CUDAcommon::getCUDAvars().gpu_lambda, sizeof(double), cudaMemcpyDeviceToHost));
std::cout<<"check "<<max1<<" "<<max2<<" "<<max3<<" "<<lambda<<" "<<gchecklambda[0]<<endl;
#if defined(CROSSCHECK) || defined(CUDAENABLED) ||defined(CUDACROSSCHECK)
        cross_checkclass::Aux=true;
#endif

        //compute new forces
        FFM.computeForces(coord, forceAux);

        //compute direction
        //TODO write the following function in CUDA
        newGrad = CGMethod::allFADotFA();
        prevGrad = CGMethod::allFADotFAP();

        //Polak-Ribieri update
        beta = max(0.0, (newGrad - prevGrad) / curGrad);

        //shift gradient
        shiftGradient(beta);

        //direction reset if not downhill, or no progress made
        if(CGMethod::allFDotFA() <= 0 || areEqual(curGrad, newGrad)) {

            shiftGradient(0.0);
            _safeMode = true;
        }
        curGrad = newGrad;

    }

    if (numIter >= N) {
        cout << endl;

        cout << "WARNING: Did not minimize in N = " << N << " steps." << endl;
//        cout << "Maximum force in system = " << maxF() << endl;

//        cout << "Culprit ..." << endl;
//        auto b = maxBead();
//        if(b != nullptr) b->getParent()->printSelf();

        cout << "System energy..." << endl;
        FFM.computeEnergy(coord, force, 0.0, true);

        cout << endl;
    }

//    cout << "Minimized." << endl;

#if defined(CROSSCHECK) || defined(CUDAENABLED) || defined(CUDACROSSCHECK)
    cross_checkclass::Aux=false;
#endif

    //final force calculation
    FFM.computeForces(coord, force);
    CUDAcommon::handleerror(cudaMemcpy(F_i, CUDAcommon::getCUDAvars().gpu_force, 3 * Bead::getBeads().size() *
                                                                                 sizeof(double), cudaMemcpyDeviceToHost));
    max2=0.0;
    std::cout<<"forces"<<endl;
    for(int iter=0;iter<Bead::getBeads().size();iter++) {
//        std::cout << F_i[3 * iter] << " " << F_i[3 * iter + 1] << " " << F_i[3 * iter + 2] <<" ";
//        std::cout <<force[3 * iter] << " " << force[3 * iter + 1] << " " << force[3 * iter + 2] << endl;
        if(abs(F_i[3 * iter]-force[3 * iter])>max2)
            max2 = abs(F_i[3 * iter]-force[3 * iter]);
        else if(abs(F_i[3 * iter +1]-force[3 * iter +1])>max2)
            max2 = abs(F_i[3 * iter +1]-force[3 * iter +1]);
        else if(abs(F_i[3 * iter +2]-force[3 * iter +2])>max2)
            max2 = abs(F_i[3 * iter +1]-force[3 * iter +1]);
    }

    std::cout<<"force alone end "<<max2<<endl;

#if defined(CROSSCHECK)
        auto state=cross_check::crosscheckforces(force);
#endif

    FFM.copyForces(forceAux, force);

    endMinimization();
    FFM.computeLoadForces();
std::cout<<"-----------------"<<endl;

    FFM.cleanupAllForceFields();
}
