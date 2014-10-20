//
//  SubSystem.cpp
//  CytoMech
//
//  Created by Konstantin Popov on 4/15/14.
//  Copyright (c) 2014 Konstantin Popov. All rights reserved.
//

#include "SubSystem.h"

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
void SubSystem::AddNewLinkers(std::vector<std::vector<Cylinder* >>& v, short linkerType){
    
    for (auto it: v) {
        //Call Linkers constructor
        LinkerDB::Instance(LinkerDBKey())->CreateLinker(it[0], it[1], linkerType);
    }
    
}
/// Add a linker
void SubSystem::AddNewLinker(Cylinder* pc1, Cylinder* pc2, short linkerType, double position1, double position2){
    
    LinkerDB::Instance(LinkerDBKey())->CreateLinker(pc1, pc2, linkerType, position1, position2);  //Call Linker constructor
}

/// Add many motor ghosts:
void SubSystem::AddNewMotorGhosts(std::vector<std::vector<Cylinder* >>& v, short motorType){
    
    for (auto it: v) {
        ///Call motors constructor
        MotorGhostDB::Instance(MotorGhostDBKey())->CreateMotorGhost(it[0], it[1], motorType);
    }

}
/// Add a motor ghost:
void SubSystem::AddNewMotorGhost(Cylinder* pc1, Cylinder* pc2, short motorType, double position1, double position2){
    
    MotorGhostDB::Instance(MotorGhostDBKey())->CreateMotorGhost(pc1, pc2, motorType, position1, position2);
    
}

// compute energy: iterate over all filaments, motors, lincers etc and call compute energy. add it up and return.
double SubSystem::getSubSystemEnergy() {return _energy;}

void SubSystem::setSubSystemEnergy(double energy) {_energy = energy;}
