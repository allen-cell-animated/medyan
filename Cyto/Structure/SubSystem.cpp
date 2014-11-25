//
//  SubSystem.cpp
//  CytoMech
//
//  Created by Konstantin Popov on 4/15/14.
//  Copyright (c) 2014 Konstantin Popov. All rights reserved.
//

#include "SubSystem.h"

#include "FilamentDB.h"
#include "LinkerDB.h"
#include "MotorGhostDB.h"


void SubSystem::addNewFilaments(vector<vector<vector<double> >>& v){
    for (auto it: v) FilamentDB::instance()->createFilament(this, it);
}
void SubSystem::addNewFilament(vector<vector<double>>& v) {
    FilamentDB::instance()->createFilament(this, v);
}
void SubSystem::removeFilament(Filament* f) {
    FilamentDB::instance()->removeFilament(f);
}

void SubSystem::addNewLinkers(vector<vector<Cylinder* >>& v, short linkerType) {
    for (auto it: v) { LinkerDB::instance()->createLinker(it[0], it[1], linkerType); }
}
void SubSystem::addNewLinker(Cylinder* c1, Cylinder* c2,
                             short linkerType, double position1, double position2){
    LinkerDB::instance()->createLinker(c1, c2, linkerType, position1, position2, true);
}
void SubSystem::removeLinker(Linker* l) {
    LinkerDB::instance()->removeLinker(l);
}

void SubSystem::addNewMotorGhosts(vector<vector<Cylinder* >>& v, short motorType) {
    for (auto it: v) MotorGhostDB::instance()->createMotorGhost(it[0], it[1], motorType);
}
void SubSystem::addNewMotorGhost(Cylinder* c1, Cylinder* c2,
                                 short motorType, double position1, double position2) {
    MotorGhostDB::instance()->createMotorGhost(c1, c2, motorType, position1, position2, true);
}
void SubSystem::removeMotorGhost(MotorGhost* m) {
    MotorGhostDB::instance()->removeMotorGhost(m);
}

double SubSystem::getSubSystemEnergy() {return _energy;}
void SubSystem::setSubSystemEnergy(double energy) {_energy = energy;}
