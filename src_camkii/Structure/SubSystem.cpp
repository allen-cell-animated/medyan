
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v3.1
//
//  Copyright (2015-2016)  Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the MEDYAN web page for more information:
//  http://www.medyan.org
//------------------------------------------------------------------

#include "SubSystem.h"

#include "BoundaryElement.h"
#include "CompartmentGrid.h"
#include "BindingManager.h"

void SubSystem::updateBindingManagers() {
    
    for(auto C : _compartmentGrid->getCompartments()) {
        
        for(auto &manager : C->getFilamentBindingManagers())
            manager->updateAllPossibleBindings();
    }
}


