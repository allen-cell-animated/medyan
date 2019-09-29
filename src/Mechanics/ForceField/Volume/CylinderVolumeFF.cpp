
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

void CylinderVolumeFF::vectorize() {

    for (auto &interaction : _cylinderVolInteractionVector)
        interaction->vectorize();
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
        /*
        cout<< "I'm here"<<endl;
        cout<<"u_i is "<< U_i<<endl;
        if(U_i > SysParams::Mechanics().cylThresh){
            cout<< "I'm here now"<<endl;
        vector<tuple<floatingpoint, floatingpoint, floatingpoint, floatingpoint, floatingpoint>> cylEnergies =  interaction->getCylEnergies();
            cout<<"size is "<<cylEnergies.size()<<endl;
            for(int i = 0; i < cylEnergies.size(); i++){

                tuple<floatingpoint, floatingpoint, floatingpoint, floatingpoint, floatingpoint> cyl = cylEnergies[i];
                cout << get<0>(cyl) << " " << get<1>(cyl) << " " << get<2>(cyl) << " " << get<3>(cyl) << " " << get<4>(cyl) <<  endl;
            };
            
            
        };
         */
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
