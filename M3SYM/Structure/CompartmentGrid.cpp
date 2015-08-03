
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

species_copy_t CompartmentGrid::countDiffusingSpecies(const string& name) {
    
    species_copy_t copyNum = 0;

    for(auto &c : children()) {
        
        auto s = ((Compartment*)(c.get()))->findSpeciesByName(name);
        assert(s != nullptr && "Counting a diffusing species that does not exist.");
        
        copyNum += s->getN();
    }
    return copyNum;
}


species_copy_t CompartmentGrid::countBulkSpecies(const string& name) {
    
    auto s = findSpeciesBulkByName(name);
    assert(s != nullptr && "Counting a bulk species that does not exist.");
    
    return s->getN();
}