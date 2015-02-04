
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

#include "Filament.h"
#include "Linker.h"
#include "MotorGhost.h"
#include "BranchingPoint.h"
#include "BoundaryElement.h"

#include "MathFunctions.h"
#include "SystemParameters.h"

//FILAMENTS
void SubSystem::addNewFilaments(vector<vector<vector<double> >>& v){
    
    for (auto it: v) {
        double d = mathfunc::twoPointDistance(it[0], it[1]);
        vector<double> tau = mathfunc::twoPointDirection(it[0], it[1]);
        
        int numSegment = d / SystemParameters::Geometry().cylinderSize;
        
        // check how many segments can fit between end-to-end of the filament
        if (numSegment == 0) new Filament(this, it[0], tau, false, false);
        else new Filament(this, it, numSegment + 1, "STRAIGHT");
    }
}
Filament* SubSystem::addNewFilament(vector<double>& position,
                                    vector<double>& direction, bool branch) {
    return new Filament(this, position, direction, true, branch);
}
void SubSystem::removeFilament(Filament* f) { delete f; }


//LINKERS
void SubSystem::addNewLinkers(vector<vector<Cylinder* >>& v,
                              short linkerType) {
    
    for (auto it: v) new Linker(it[0], it[1], linkerType);
}
Linker* SubSystem::addNewLinker(Cylinder* c1, Cylinder* c2,
                                short linkerType,
                                double position1,
                                double position2){
    return new Linker(c1, c2, linkerType, position1, position2, true);
}
void SubSystem::removeLinker(Linker* l) {delete l;}


//MOTORGHOSTS
void SubSystem::addNewMotorGhosts(vector<vector<Cylinder* >>& v,
                                  short motorType) {
    
    for (auto it: v) new MotorGhost(it[0], it[1], motorType);
}
MotorGhost* SubSystem::addNewMotorGhost(Cylinder* c1, Cylinder* c2,
                                        short motorType,
                                        double position1,
                                        double position2) {
    return new MotorGhost(c1, c2, motorType, position1, position2, true);
}
void SubSystem::removeMotorGhost(MotorGhost* m) { delete m; }



//BRANCHINGPOINTS
void SubSystem::addNewBranchingPoints(vector<vector<Cylinder* >>& v,
                                      short branchType) {
    
    for (auto it: v) new BranchingPoint(it[0], it[1], branchType);
}
BranchingPoint* SubSystem::addNewBranchingPoint(Cylinder* c1, Cylinder* c2,
                                                short branchType,
                                                double position) {
    return new BranchingPoint(c1, c2, branchType, position, true);
}
void SubSystem::removeBranchingPoint(BranchingPoint* b) { delete b; }


//OTHER SUBSYSTEM FUNCTIONS
double SubSystem::getSubSystemEnergy() {return _energy;}
void SubSystem::setSubSystemEnergy(double energy) {_energy = energy;}

#ifdef DYNAMICRATES
vector<Cylinder*> SubSystem::getBoundaryCylinders() {

    vector<Cylinder*> cylinders;
    auto list = getNeighborList();
    
    //loop through neighbor list, construct vector
    for(auto be : *BoundaryElementDB::instance()) {
        auto localCylinders = list->getNeighbors(be);
        cylinders.insert(cylinders.end(), localCylinders.begin(),
                                          localCylinders.end());
    }
    return cylinders;
}
#endif


