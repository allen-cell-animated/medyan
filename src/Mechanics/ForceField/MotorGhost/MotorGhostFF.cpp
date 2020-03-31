
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

#include "MotorGhostFF.h"

#include "MotorGhostStretching.h"
#include "MotorGhostStretchingHarmonic.h"

#include "MotorGhost.h"

void MotorGhostFF::assignforcemags(){
    for (auto &interaction : _motorGhostInteractionVector)
        interaction->assignforcemags();
}
MotorGhostFF::MotorGhostFF (string& stretching, string& bending, string& twisting)
{
    if (stretching == "HARMONIC")
        _motorGhostInteractionVector.emplace_back(
                new MotorGhostStretching<MotorGhostStretchingHarmonic>());
    else if(stretching == "") {}
    else {
        cout << "Motor stretching FF not recognized. Exiting." << endl;
        exit(EXIT_FAILURE);
    }
}

void MotorGhostFF::whoIsCulprit() {

    cout << endl;

    cout << "Culprit interaction = " << _culpritInteraction->getName() << endl;

    cout << "Printing the culprit motor..." << endl;
    _culpritInteraction->_motorCulprit->printSelf();

    cout << endl;
}

void MotorGhostFF::vectorize(const FFCoordinateStartingIndex& si) {
    //Reset stretching forces to 0.
    for(auto m:MotorGhost::getMotorGhosts()){
        //Using += to ensure that the stretching forces are additive.
        m->getMMotorGhost()->stretchForce = 0.0;
    }


    for (auto &interaction : _motorGhostInteractionVector)
        interaction->vectorize(si);
}

void MotorGhostFF::cleanup() {

    for (auto &interaction : _motorGhostInteractionVector)
        interaction->deallocate();
}


floatingpoint MotorGhostFF::computeEnergy(floatingpoint *coord, bool stretched) {

    floatingpoint U= 0.0;
    floatingpoint U_i=0.0;

    for (auto &interaction : _motorGhostInteractionVector) {

        U_i = interaction->computeEnergy(coord);

        if(U_i <= -1) {
            //set culprit and return
            _culpritInteraction = interaction.get();
            return -1;
        }
        else U += U_i;

        #ifdef TRACKDIDNOTMINIMIZE
        if(!stretched)
            SysParams::Mininimization().tempEnergyvec.push_back(U_i);
        #endif
    }
    return U;
}

void MotorGhostFF::computeForces(floatingpoint *coord, floatingpoint *f) {

    for (auto &interaction : _motorGhostInteractionVector)
        interaction->computeForces(coord, f);
}

vector<string> MotorGhostFF::getinteractionnames(){
	vector<string> temp;
	for (auto &interaction : _motorGhostInteractionVector) {
		temp.push_back(interaction->getName());
	}
    return temp;
}


