
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

#include "CompartmentGrid.h"

#include "ChemSim.h"

void CompartmentGrid::addChemSimReactions(ChemSim* chem) {
    
    for(auto &c : children()) {
        Compartment* C = (Compartment*)(c.get());
        C->addChemSimReactions(chem);
    }
    
    for(auto &r : _bulkReactions.reactions())
        chem->addReaction(r.get());
    
}