
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

#include "CylinderVolumeFF.h"

#include "CylinderExclVolume.h"
#include "CylinderExclVolRepulsion.h"

#include "Cylinder.h"

CylinderVolumeFF::CylinderVolumeFF (string& type) {
    if (type == "REPULSION")
        _cylinderVolInteractionVector.emplace_back(
        new CylinderExclVolume <CylinderExclVolRepulsion>());
    else if(type == "") {}
    else {
        cout << "Volume FF not recognized. Exiting." << endl;
        exit(EXIT_FAILURE);
    }
}

void CylinderVolumeFF::whoIsCulprit() {

    cout << endl;

    cout << "Culprit interaction = " << _culpritInteraction->getName() << endl;

    cout << "Printing the culprit cylinders..." << endl;
    _culpritInteraction->_cylinderCulprit1->printSelf();
    _culpritInteraction->_cylinderCulprit2->printSelf();

    cout << endl;
}

void CylinderVolumeFF::vectorize(const FFCoordinateStartingIndex& si) {

    for (auto &interaction : _cylinderVolInteractionVector)
        interaction->vectorize(si);
}

void CylinderVolumeFF::cleanup() {

    for (auto &interaction : _cylinderVolInteractionVector)
        interaction->deallocate();
}


floatingpoint CylinderVolumeFF::computeEnergy(floatingpoint *coord, bool stretched) {

    floatingpoint U= 0.0;
    floatingpoint U_i=0.0;

    for (auto &interaction : _cylinderVolInteractionVector) {

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

#if defined(HYBRID_NLSTENCILLIST) || defined(SIMDBINDINGSEARCH)
void CylinderVolumeFF::setHNeighborLists(HybridCylinderCylinderNL* Hnl) {
    for (auto &interaction : _cylinderVolInteractionVector){
        interaction->setHNeighborList(Hnl);
    }
};
#endif

void CylinderVolumeFF::computeForces(floatingpoint *coord, floatingpoint *f) {

    for (auto &interaction : _cylinderVolInteractionVector)
        interaction->computeForces(coord, f);
}


vector<NeighborList*> CylinderVolumeFF::getNeighborLists() {

    vector<NeighborList*> neighborLists;

    for(auto &interaction : _cylinderVolInteractionVector)
        neighborLists.push_back(interaction->getNeighborList());

    return neighborLists;
}

vector<string> CylinderVolumeFF::getinteractionnames(){
	vector<string> temp;
	for (auto &interaction : _cylinderVolInteractionVector) {
		temp.push_back(interaction->getName());
	}
    return temp;
}

