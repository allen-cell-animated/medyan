
//------------------------------------------------------------------
//  **M3SYM** - Simulation Package for the Mechanochemical
//              Dynamics of Active Networks, 3rd Generation
//
//  Copyright (2014) Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the Papoian lab page for more information:
//  http://papoian.chem.umd.edu/
//------------------------------------------------------------------

#include "BoundaryFF.h"

#include "BoundaryRepulsion.h"
#include "BoundaryRepulsionLJ.h"
#include "BoundaryRepulsionExp.h"

#include "BoundaryElement.h"
#include "Bead.h"

BoundaryFF::BoundaryFF (string type) {
    
    if (type == "REPULSIONLJ")
        _boundaryInteractionVector.emplace_back(
            new BoundaryRepulsion<BoundaryRepulsionLJ>());
    else if (type == "REPULSIONEXP")
        _boundaryInteractionVector.emplace_back(
            new BoundaryRepulsion<BoundaryRepulsionExp>());
    else if(type == "") {}
    else {
        cout << "Boundary FF not recognized. Exiting." << endl;
        exit(EXIT_FAILURE);
    }
}

void BoundaryFF::whoIsCulprit() {
    
    cout << endl;
    
    cout << "Printing the culprit bead and boundary element..." << endl;
    
    _beadCulprit->printInfo();
    _boundaryCulprit->printInfo();
    
    cout << endl;
}


double BoundaryFF::computeEnergy(double d) {
    
    double U = 0;
    double U_i;
    
    for (auto &interaction : _boundaryInteractionVector){
        
        auto nl = interaction->getNeighborList();
        
        for (auto be: BoundaryElement::getBoundaryElements()) {
            
            for(auto &bd : nl->getNeighbors(be)) {
                
                U_i = interaction->computeEnergy(be, bd, d);
                
                if(fabs(U_i) == numeric_limits<double>::infinity() || U_i != U_i) {
                    
                    //set culprits and return
                    _beadCulprit = bd;
                    _boundaryCulprit = be;
                    
                    return -1;
                }
                else
                    U += U_i;
            }
        }
    }
    return U;
}

void BoundaryFF::computeForces() {

    for (auto &interaction : _boundaryInteractionVector){
        
        auto nl = interaction->getNeighborList();
        
        for (auto be: BoundaryElement::getBoundaryElements()) {
            
            for(auto bd : nl->getNeighbors(be))
                interaction->computeForces(be, bd);
        }
    }
}

void BoundaryFF::computeForcesAux() {
    
    for (auto &interaction : _boundaryInteractionVector){
        
        auto nl = interaction->getNeighborList();
        
        for (auto be: BoundaryElement::getBoundaryElements()) {
            
            for(auto bd : nl->getNeighbors(be))
                interaction->computeForcesAux(be, bd);
        }
    }
}

vector<NeighborList*> BoundaryFF::getNeighborLists() {
    
    vector<NeighborList*> neighborLists;
    
    for(auto &interaction : _boundaryInteractionVector)
        neighborLists.push_back(interaction->getNeighborList());
    
    return neighborLists;
}
