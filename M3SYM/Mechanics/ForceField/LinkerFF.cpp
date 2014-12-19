
//------------------------------------------------------------------
//  **M3SYM** - Simulation Package for the Mechanochemical
//              Dynamics of Active Networks, 3rd Generation
//
//  Copyright (2014) Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the Papoian lab page for installation and documentation:
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
}

double LinkerFF::computeEnergy(double d) {
    double U_linker = 0;
    
    for (auto linker: *LinkerDB::instance())
        for (auto &linkerInteraction : _linkerInteractionVector)
            U_linker += linkerInteraction.get()->computeEnergy(linker, d);

    return U_linker;
}

void LinkerFF::computeForces() {
    
    for (auto linker: *LinkerDB::instance())
        for (auto &linkerInteraction : _linkerInteractionVector)
            linkerInteraction.get()->computeForces(linker);
}

void LinkerFF::computeForcesAux() {
    
    for (auto linker: *LinkerDB::instance())
        for (auto &linkerInteraction : _linkerInteractionVector)
            linkerInteraction.get()->computeForcesAux(linker);
}

