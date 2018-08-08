
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

#include "FilamentFF.h"

#include "FilamentStretching.h"
#include "FilamentStretchingHarmonic.h"

#include "FilamentBending.h"
#include "FilamentBendingHarmonic.h"
#include "FilamentBendingCosine.h"

#include "Filament.h"
//TODO remove later
#include "CGMethod.h"
#include "cross_check.h"

FilamentFF::FilamentFF (string& stretching, string& bending, string& twisting) {

    if (stretching == "HARMONIC")
        _filamentInteractionVector.emplace_back(
                new FilamentStretching<FilamentStretchingHarmonic>());
    else if(stretching == "") {}
    else {
        cout << "Filament stretching FF not recognized. Exiting." << endl;
        exit(EXIT_FAILURE);
    }

    if (bending == "HARMONIC")
        _filamentInteractionVector.emplace_back(
                new FilamentBending<FilamentBendingHarmonic>());
    else if(bending == "COSINE")
        _filamentInteractionVector.emplace_back(
                new FilamentBending<FilamentBendingCosine>());
    else if(bending == "") {}
    else {
        cout << "Filament bending FF not recognized. Exiting." << endl;
        exit(EXIT_FAILURE);
    }
}

void FilamentFF::vectorize() {

    for (auto &interaction : _filamentInteractionVector)
        interaction->vectorize();
}

void FilamentFF::cleanup() {

    for (auto &interaction : _filamentInteractionVector)
        interaction->deallocate();
}

void FilamentFF::whoIsCulprit() {

    cout << endl;

    cout << "Culprit interaction = " << _culpritInteraction->getName() << endl;

    cout << "Printing the culprit filament..." << endl;
    _culpritInteraction->_filamentCulprit->printSelf();

    cout << endl;
}


double FilamentFF::computeEnergy(double *coord, double *f, double d) {

    double U= 0.0;
    double U_i=0.0;

    for (auto &interaction : _filamentInteractionVector) {
//        std::cout<<"ForceField "<<interaction->getName()<<" "<<_filamentInteractionVector.size()<<endl;
        U_i = interaction->computeEnergy(coord, f, d);
//        CUDAcommon::handleerror(cudaDeviceSynchronize(),"FilamentFF","FilamentFF");
        if(U_i <= -1) {
            //set culprit and return
            _culpritInteraction = interaction.get();
            return -1;
        }
        else U += U_i;
        
        if(U_i > 1000.0){
            std::cout<<" High Energy in" <<getName()<<" = "<<U_i<<endl;
        }
#ifdef DETAILEDOUTPUT
        std::cout<<getName()<<" "<<U_i<<endl;
#endif
    }

    return U;
}

void FilamentFF::computeForces(double *coord, double *f) {
//    double *F_i = new double[CGMethod::N];
    for (auto &interaction : _filamentInteractionVector) {
//        std::cout<<"ForceField "<<interaction->getName()<<" "<<_filamentInteractionVector.size()<<endl;
        interaction->computeForces(coord, f);
//        CUDAcommon::handleerror(cudaDeviceSynchronize());


//        if(cross_checkclass::Aux)
//            CUDAcommon::handleerror(
//                    cudaMemcpy(F_i, CUDAcommon::getCUDAvars().gpu_forceAux, CGMethod::N * sizeof
//                                       (double),
//                               cudaMemcpyDeviceToHost));
//        else
//            CUDAcommon::handleerror(
//                    cudaMemcpy(F_i, CUDAcommon::getCUDAvars().gpu_force, CGMethod::N * sizeof
//                                       (double),
//                               cudaMemcpyDeviceToHost));
//        double fmax = 0.0;
//        int id=0;
//        for (auto iter = 0; iter < CGMethod::N/3; iter++) {
//            if(abs(F_i[3 *iter])> fmax) {fmax = abs(F_i[3*iter]);id = iter;}
//            if(abs(F_i[3 *iter +1])> fmax) {fmax = abs(F_i[3*iter +1]);id = iter;}
//            if(abs(F_i[3 *iter +2])> fmax) {fmax = abs(F_i[3*iter +2]);id = iter;}
////            std::cout << F_i[3 * iter] << " " << F_i[3 * iter + 1] << " " << F_i[3 * iter + 2] << endl;
//        }
//        std::cout <<"Fmax "<< id<<" "<<fmax<<" "<<F_i[3 * id] << " " << F_i[3 * id + 1] << " " << F_i[3 * id + 2] <<
//                  endl;
    }
    //TODO remove later
//    delete F_i;
}
