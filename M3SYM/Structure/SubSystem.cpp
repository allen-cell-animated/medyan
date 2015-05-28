
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

#include "SubSystem.h"

#include "BoundaryElement.h"
#include "CompartmentGrid.h"
#include "BindingManager.h"

void SubSystem::updateBindingManagers() {
    for(auto &child : _compartmentGrid->children()) {
        Compartment* c = (Compartment*)child.get();
        
        for(auto &manager : c->getFilamentBindingManagers())
            manager->updateAllPossibleBindings();
    }
}


