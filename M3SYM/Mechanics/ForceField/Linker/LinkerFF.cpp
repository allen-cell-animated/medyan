
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

#include "LinkerFF.h"

#include "LinkerStretching.h"
#include "LinkerStretchingHarmonic.h"

#include "Linker.h"

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

double LinkerFF::computeEnergy(double d) {
    
    double U = 0;
    double U_i;
    
    for (auto &interaction : _linkerInteractionVector) {
        
        for (auto l: Linker::getLinkers()) {
            
            U_i = interaction->computeEnergy(l, d);
            
            if(fabs(U_i) == numeric_limits<double>::infinity() || U_i != U_i)
                return -1;
            else
                U += U_i;
        }
    }
    return U;
}

void LinkerFF::computeForces() {
    
    for (auto &interaction : _linkerInteractionVector) {
        
        for (auto l: Linker::getLinkers())
            interaction->computeForces(l);
    }
}

void LinkerFF::computeForcesAux() {
    
    for (auto &interaction : _linkerInteractionVector) {
        
        for (auto l: Linker::getLinkers())
            interaction->computeForcesAux(l);
    }
}

