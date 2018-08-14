
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v3.2
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

//vector<tuple<int, int, vector<double>, vector<double>>> CompartmentGrid::getDiffusingMotors() {
//    
//    vector<tuple<int, int, vector<double>, vector<double>>> output;
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
//                    double dist = (mm->getRMax() + mm->getRMin()) / 2.0;
//                
//                    vector<double> midpoint = GController::getRandomCoordinates(c);
//                
//                    double directionX = Rand::randDouble(-1,1);
//                    double directionY = Rand::randDouble(-1,1);
//                    double directionZ = Rand::randDouble(-1,1);
//                    vector<double> direction = normalizedVector({directionX, directionY, directionZ});
//                
//                    vector<double> firstPoint  = nextPointProjection(midpoint, dist / 2.0, direction);
//                    vector<double> secondPoint = nextPointProjection(midpoint, dist / 2.0,
//                                                 vector<double>{-direction[0], -direction[1], -direction[2]});
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

