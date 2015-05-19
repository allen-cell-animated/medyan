
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
        _BoundaryInteractionVector.emplace_back(
            new BoundaryRepulsion<BoundaryRepulsionLJ>());
    else if (type == "REPULSIONEXP")
        _BoundaryInteractionVector.emplace_back(
            new BoundaryRepulsion<BoundaryRepulsionExp>());
    else if(type == "") {}
    else {
        cout << "Boundary FF not recognized. Exiting." << endl;
        exit(EXIT_FAILURE);
    }
}

double BoundaryFF::computeEnergy(double d) {
    
    double U = 0;
    double U_i;
    
    for (auto &interaction : _BoundaryInteractionVector){
        
        auto nl = interaction->getNeighborList();
        
        for (auto be: BoundaryElement::getBoundaryElements()) {
            
            for(auto &bd : nl->getNeighbors(be)) {
                
                U_i = interaction->computeEnergy(be, bd, d);
                
                if(fabs(U_i) == numeric_limits<double>::infinity() || U_i != U_i)
                    return -1;
                else
                    U += U_i;
            }
        }
    }
    return U;
}

void BoundaryFF::computeForces() {

    for (auto &interaction : _BoundaryInteractionVector){
        
        auto nl = interaction->getNeighborList();
        
        for (auto be: BoundaryElement::getBoundaryElements()) {
            
            for(auto bd : nl->getNeighbors(be))
                interaction->computeForces(be, bd);
        }
    }
}

void BoundaryFF::computeForcesAux() {
    
    for (auto &interaction : _BoundaryInteractionVector){
        
        auto nl = interaction->getNeighborList();
        
        for (auto be: BoundaryElement::getBoundaryElements()) {
            
            for(auto bd : nl->getNeighbors(be))
                interaction->computeForcesAux(be, bd);
        }
    }
}
