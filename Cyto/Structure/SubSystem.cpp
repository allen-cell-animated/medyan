//
//  SubSystem.cpp
//  CytoMech
//
//  Created by Konstantin Popov on 4/15/14.
//  Copyright (c) 2014 Konstantin Popov. All rights reserved.
//

#include "SubSystem.h"

#include "BeadDB.h"
#include "CylinderDB.h"
#include "FilamentDB.h"
#include "MotorGhostDB.h"
#include "LinkerDB.h"

using namespace std;

// Interfaces for creation of a Filament (vector with coordinates of beginning and end for the filament)
void SubSystem::AddNewFilaments(vector<vector<vector<double> >>& v){
    
    for (auto it: v) FilamentDB::Instance(FilamentDBKey())->CreateFilament(this, it);
    //Call FilamentDB constructor
}

/// Add many linkers:
void SubSystem::AddNewLinkers(std::vector<std::vector<Cylinder* >>& v){
    
    for (auto it: v) {
        ///find compartment
        auto m1 = MidPointCoordinate(it[0]->GetFirstBead()->coordinate, it[0]->GetSecondBead()->coordinate, 0.5);
        auto m2 = MidPointCoordinate(it[1]->GetFirstBead()->coordinate, it[1]->GetSecondBead()->coordinate, 0.5);
        auto position = MidPointCoordinate(m1, m2, 0.5);
        
        Compartment* c;
        try {c = GController::getCompartment(position);}
        catch (std::exception& e) {std:: cout << e.what(); exit(EXIT_FAILURE);}
        
        //Call Linkers constructor
        LinkerDB::Instance(LinkerDBKey())->CreateLinker(it[0], it[1], c);
    }
    
}
/// Add a linker
void SubSystem::AddNewLinker(Cylinder* pc1, Cylinder* pc2, Compartment* c, double position1, double position2){
    
    LinkerDB::Instance(LinkerDBKey())->CreateLinker(pc1, pc2, c, position1, position2);  //Call Linker constructor
}

/// Add many motor ghosts:
void SubSystem::AddNewMotorGhosts(std::vector<std::vector<Cylinder* >>& v, double k, double position1, double position2){
    
    for (auto it: v)
        MotorGhostDB::Instance(MotorGhostDBKey())->CreateMotorGhost(it[0], it[1], k, position1, position2);

}
/// Add a motor ghost:
void SubSystem::AddNewMotorGhost(Cylinder* pc1, Cylinder* pc2, double k, double position1, double position2){
    
    MotorGhostDB::Instance(MotorGhostDBKey())->CreateMotorGhost(pc1, pc2, k, position1, position2);
    
}

// compute energy: iterate over all filaments, motors, lincers etc and call compute energy. add it up and return.
double SubSystem::getSubSystemEnergy() {return _energy;}

void SubSystem::setSubSystemEnergy(double energy) {_energy = energy;}
