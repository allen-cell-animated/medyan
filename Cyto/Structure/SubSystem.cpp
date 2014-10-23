//
//  SubSystem.cpp
//  CytoMech
//
//  Created by Konstantin Popov on 4/15/14.
//  Copyright (c) 2014 Konstantin Popov. All rights reserved.
//

#include "SubSystem.h"

using namespace std;

void SubSystem::AddNewFilaments(vector<vector<vector<double> >>& v){
    
    for (auto it: v) FilamentDB::Instance(FilamentDBKey())->CreateFilament(this, it);
}

void SubSystem::AddNewFilament(vector<vector<double>>& v) {
    
    FilamentDB::Instance(FilamentDBKey())->CreateFilament(this, v);
}

void SubSystem::RemoveFilament(Filament* f) {
    
    FilamentDB::Instance(FilamentDBKey())->RemoveFilament(f);
}

void SubSystem::AddNewLinkers(std::vector<std::vector<Cylinder* >>& v, short linkerType){
    
    for (auto it: v) { LinkerDB::Instance(LinkerDBKey())->CreateLinker(it[0], it[1], linkerType); }
    
}
void SubSystem::AddNewLinker(Cylinder* pc1, Cylinder* pc2, short linkerType, double position1, double position2){
    
    LinkerDB::Instance(LinkerDBKey())->CreateLinker(pc1, pc2, linkerType, position1, position2, true);
}

void SubSystem::RemoveLinker(Linker* l) {
    
    LinkerDB::Instance(LinkerDBKey())->RemoveLinker(l);
}

void SubSystem::AddNewMotorGhosts(std::vector<std::vector<Cylinder* >>& v, short motorType){
    
    for (auto it: v) { MotorGhostDB::Instance(MotorGhostDBKey())->CreateMotorGhost(it[0], it[1], motorType); }

}
void SubSystem::AddNewMotorGhost(Cylinder* pc1, Cylinder* pc2, short motorType, double position1, double position2){
    
    MotorGhostDB::Instance(MotorGhostDBKey())->CreateMotorGhost(pc1, pc2, motorType, position1, position2, true);
}

void SubSystem::RemoveMotorGhost(MotorGhost* m) {
    
    MotorGhostDB::Instance(MotorGhostDBKey())->RemoveMotorGhost(m);
}

double SubSystem::getSubSystemEnergy() {return _energy;}

void SubSystem::setSubSystemEnergy(double energy) {_energy = energy;}
