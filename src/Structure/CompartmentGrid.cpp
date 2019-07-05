
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

#include "CompartmentGrid.h"
#include "ChemSim.h"

#include "MathFunctions.h"
#include "SysParams.h"
#include "GController.h"

using namespace mathfunc;

void CompartmentGrid::addChemSimReactions(ChemSim* chem) {
    
    for(auto C : getCompartments())
        C->addChemSimReactions(chem);
    
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

//DEPRECATED AS OF 9/8/16

//vector<tuple<int, int, vector<floatingpoint>, vector<floatingpoint>>> CompartmentGrid::getDiffusingMotors() {
//    
//    vector<tuple<int, int, vector<floatingpoint>, vector<floatingpoint>>> output;
//
//    //for all motor types, get compartment binding managers
//    //@note - this is for filament type 0 only. For multiple filament types, this would require a fix.
//    for(int type = 0; type < SysParams::Chemistry().numMotorSpecies[0]; type++) {
//        
//        for(Compartment* c: getCompartments()) {
//            
//            MotorBindingManager* mm = c->getMotorBindingManager(type);
//            
//            for(int ID : mm->getAllUnboundIDs()) {
//            
//                bool found = false;
//                
//                while(!found) {
//                
//                    //pick random two points in compartment
//                    floatingpoint dist = (mm->getRMax() + mm->getRMin()) / 2.0;
//                
//                    vector<floatingpoint> midpoint = GController::getRandomCoordinates(c);
//                
//                    floatingpoint directionX = Rand::randfloatingpoint(-1,1);
//                    floatingpoint directionY = Rand::randfloatingpoint(-1,1);
//                    floatingpoint directionZ = Rand::randfloatingpoint(-1,1);
//                    vector<floatingpoint> direction = normalizeVector({directionX, directionY, directionZ});
//                
//                    vector<floatingpoint> firstPoint  = nextPointProjection(midpoint, dist / 2.0, direction);
//                    vector<floatingpoint> secondPoint = nextPointProjection(midpoint, dist / 2.0,
//                                                 vector<floatingpoint>{-direction[0], -direction[1], -direction[2]});
//                    
//                    try {
//                        GController::getCompartment(firstPoint);
//                        GController::getCompartment(secondPoint);
//                        
//                        //add to output and switch flag
//                        output.emplace_back(ID, type, firstPoint, secondPoint);
//                        found = true;
//                    }
//                    catch (exception& e) {/*just repeat loop*/}
//                }
//                
//            }
//        }
//    }
//    return output;
//}

