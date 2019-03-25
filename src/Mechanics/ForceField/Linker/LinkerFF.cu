
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

#include "LinkerFF.h"

#include "LinkerStretching.h"
#include "LinkerStretchingHarmonic.h"

#include "Linker.h"

void LinkerFF::assignforcemags(){
    for (auto &interaction : _linkerInteractionVector)
        interaction->assignforcemags();
}

LinkerFF::LinkerFF (string& stretching, string& bending, string& twisting)
{
    if (stretching == "HARMONIC")
        _linkerInteractionVector.emplace_back(
                new LinkerStretching<LinkerStretchingHarmonic>());
    else if(stretching == "") {}
    else {
        cout << "Linker stretching FF not recognized. Exiting." << endl;
        exit(EXIT_FAILURE);
    }
}

void LinkerFF::vectorize() {
    //Reset stretching forces to 0.
    for(auto l:Linker::getLinkers()){
        //Using += to ensure that the stretching forces are additive.
        l->getMLinker()->stretchForce = 0.0;
    }

    for (auto &interaction : _linkerInteractionVector)
        interaction->vectorize();
}

void LinkerFF::cleanup() {

    for (auto &interaction : _linkerInteractionVector)
        interaction->deallocate();
}

void LinkerFF::whoIsCulprit() {

    cout << endl;

    cout << "Culprit interaction = " << _culpritInteraction->getName() << endl;

    cout << "Printing the culprit linker..." << endl;
    _culpritInteraction->_linkerCulprit->printSelf();

    cout << endl;
}

totalenergyfloatingpoint LinkerFF::computeEnergy(floatingpoint *coord, totalforcefloatingpoint *f, floatingpoint d) {

    totalenergyfloatingpoint U= 0.0;
    floatingpoint U_i=0.0;

    for (auto &interaction : _linkerInteractionVector) {

        U_i = interaction->computeEnergy(coord, f, d);

        if(U_i <= -1) {
            //set culprit and return
            _culpritInteraction = interaction.get();
            return -1;
        }
        else U += U_i;

    }
    return U;
}

void LinkerFF::computeForces(floatingpoint *coord, totalforcefloatingpoint *f) {

    for (auto &interaction : _linkerInteractionVector)
        interaction->computeForces(coord, f);
}


