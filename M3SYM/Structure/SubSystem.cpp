
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

#include "SubSystem.h"

#include "Filament.h"
#include "Linker.h"
#include "MotorGhost.h"
#include "BranchingPoint.h"

#include "MathFunctions.h"
#include "SystemParameters.h"

void SubSystem::addNewFilaments(vector<vector<vector<double> >>& v){
    
    for (auto it: v) {
        double d = mathfunc::TwoPointDistance(it[0], it[1]);
        vector<double> tau = mathfunc::TwoPointDirection(it[0], it[1]);
        
        int numSegment = d / SystemParameters::Geometry().cylinderSize;
        
        // check how many segments can fit between end-to-end of the filament
        if (numSegment == 0) new Filament(this, it[0], tau);
        else new Filament(this, it, numSegment + 1, "STRAIGHT");
    }
}

void SubSystem::addNewFilament(vector<vector<double>>& v) {
    
    double d = mathfunc::TwoPointDistance(v[0], v[1]);
    vector<double> tau = mathfunc::TwoPointDirection(v[0], v[1]);
    
    int numSegment = d / SystemParameters::Geometry().cylinderSize;
    
    // check how many segments can fit between end-to-end of the filament
    if (numSegment == 0) new Filament(this, v[0], tau);
    else new Filament(this, v, numSegment + 1, "STRAIGHT");
}
void SubSystem::removeFilament(Filament* f) { delete f; }


void SubSystem::addNewLinkers(vector<vector<Cylinder* >>& v, short linkerType) {
    
    for (auto it: v)
        new Linker(it[0], it[1], linkerType);
}
void SubSystem::addNewLinker(Cylinder* c1, Cylinder* c2,
                             short linkerType, double position1, double position2){
    
    new Linker(c1, c2, linkerType, position1, position2, true);
}
void SubSystem::removeLinker(Linker* l) {delete l;}


void SubSystem::addNewMotorGhosts(vector<vector<Cylinder* >>& v, short motorType) {
    
    for (auto it: v)
        new MotorGhost(it[0], it[1], motorType);
}
void SubSystem::addNewMotorGhost(Cylinder* c1, Cylinder* c2,
                                 short motorType, double position1, double position2) {
    
    new MotorGhost(c1, c2, motorType, position1, position2, true);
}
void SubSystem::removeMotorGhost(MotorGhost* m) { delete m; }

void SubSystem::addNewBranchingPoints(vector<vector<Cylinder* >>& v, short branchType) {
    
    for (auto it: v)
        new BranchingPoint(it[0], it[1], branchType);
}
void SubSystem::addNewBranchingPoint(Cylinder* c1, Cylinder* c2,
                                 short branchType, double position) {
    
    new BranchingPoint(c1, c2, branchType, position, true);
}
void SubSystem::removeBranchingPoint(BranchingPoint* b) { delete b; }


double SubSystem::getSubSystemEnergy() {return _energy;}
void SubSystem::setSubSystemEnergy(double energy) {_energy = energy;}

