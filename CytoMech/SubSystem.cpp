//
//  SubSystem.cpp
//  CytoMech
//
//  Created by Konstantin Popov on 4/15/14.
//  Copyright (c) 2014 Konstantin Popov. All rights reserved.
//

#include "SubSystem.h"

using namespace std;

//void SubSystem::AddNewFilament(vector<double> v) { _pfilamentVector.push_back( _pfdb->CreateFilament(this, v));}


// Interfaces for creation of a Filament (vector with coordinates of beginning and end for the filament)
void SubSystem::AddNewFilaments(vector<vector<vector<double> > > v){
    
    for (auto it: v) FilamentDB::Instance(FilamentDBKey())->CreateFilament(this, it);
    //Call FilamentDB constructor
}

/// Add many linkers:
void SubSystem::AddNewLinkers(std::vector<std::vector<Cylinder* > > v, double stretchConst){
    
    for (auto it: v) LinkerDB::Instance(LinkerDBKey())->CreateLinker(this, it[0], it[1], stretchConst);
    //Call Linkers constructor
    
}
/// Add a linker
void SubSystem::AddNewLinker(Cylinder* pc1, Cylinder* pc2, double stretchConst){
    
    LinkerDB::Instance(LinkerDBKey())->CreateLinker(this, pc1, pc2, stretchConst);  //Call Linker constructor
}

/// Add many motor ghosts:
void SubSystem::AddNewMotorGhosts(std::vector<std::vector<Cylinder* > > v, double k, double position1, double position2){
    
    for (auto it: v)
        MotorGhostDB::Instance(MotorGhostDBKey())->CreateMotorGhost(this, it[0], it[1], k, position1, position2);

}
/// Add a motor ghost:
void SubSystem::AddNewMotorGhost(Cylinder* pc1, Cylinder* pc2, double k, double position1, double position2){
    
    MotorGhostDB::Instance(MotorGhostDBKey())->CreateMotorGhost(this, pc1, pc2, k, position1, position2);
    
}

// compute energy: iterate over all filaments, motors, lincers etc and call compute energy. add it up and return.
double SubSystem::getSubSystemEnergy() {return _energy;}
