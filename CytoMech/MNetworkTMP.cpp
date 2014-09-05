//
//  MNetwork.cpp
//  CytoMech
//
//  Created by Konstantin Popov on 4/15/14.
//  Copyright (c) 2014 Konstantin Popov. All rights reserved.
//

#include "MNetwork.h"
#include "MSystem.h"

using namespace std;

Network::Network(System* ps)
{
    _pS = ps;
    _pfdb = new FilamentDB();
    _pldb = new LinkerDB();
    _pmgdb = new MotorGhostDB();
    
    
	
}

void Network::AddNewFilaments(vector<vector<vector<double> > > v)
{
    
    for (auto it: v){
        
        _pfdb->CreateFilament(_pS, this, it);  //Call Filament constructor
    }
}
/// Add many linkers:
void Network::AddNewLinkers(std::vector<std::vector<Cylinder* > > v, double stretchConst){
    for (auto it: v){
        _pldb->CreateLinker(this, it[0], it[1], stretchConst); //Call Linkers constructor
    }
}
/// Add a linker
void Network::AddNewLinker(Cylinder* pc1, Cylinder* pc2, double stretchConst){
    
    _pldb->CreateLinker(this, pc1, pc2, stretchConst);  //Call Linker constructor
}
/// Add many linkers:
void Network::AddNewMotorGhosts(std::vector<std::vector<Cylinder* > > v, double k, double position1, double position2){
    for (auto it: v){
        _pmgdb->CreateMotorGhost(this, it[0], it[1], k, position1, position2);
    }
}
/// Add a gost motor:
void Network::AddNewMotorGhost(Cylinder* pc1, Cylinder* pc2, double k, double position1, double position2){
    
    _pmgdb->CreateMotorGhost(this, pc1, pc2, k, position1, position2);
    
}

// compute energy: iterate over all filaments, motors, lincers etc and call compute energy. add it up ant return.
double Network::CopmuteEnergy(double d)
{
    double totEnergy = 0;
    
    if (d == 0.0){
        
        double filamentEnergy = 0;
        double motorEnergy = 0;
        double linkerEnergy = 0;
        
        for (auto it: (*_pfdb) ) {filamentEnergy += it->CopmuteEnergy();}
        
        for (auto it: (*_pmgdb) ) {motorEnergy += it->CopmuteEnergy();}
        
        for (auto it: (*_pldb) ) {linkerEnergy += it->CopmuteEnergy();}
        
        return totEnergy = filamentEnergy + linkerEnergy + motorEnergy;
    }
    
    
    else {
        
        double filamentEnergy = 0;
        double motorEnergy = 0;
        double linkerEnergy = 0;
        
        for (auto it: (*_pfdb) ) {filamentEnergy += it->CopmuteEnergy(d);}
        
        for (auto it: (*_pmgdb) ) {motorEnergy += it->CopmuteEnergy(d);}
        
        for (auto it: (*_pldb) ) {linkerEnergy += it->CopmuteEnergy(d);}
        
        return totEnergy = filamentEnergy;
    }
}

// End Energy calculation

// Force Calculation:


void Network::ResetForces(BeadDB&  list){
    for(auto it: list) {
        it->force.assign (3, 0); //Set force to zero;
    }
    
}

void Network::ResetForcesAux(BeadDB&  list){
    for(auto it: list) {
        it->forceAux.assign (3, 0); //Set forceAux to zero;
        
    }
    
}


void Network::CopmuteForce(int i)
{
    if (i == 0){
        ResetForces(*_pS->getBDB());
        for (auto it: (*_pfdb) ) {it->CopmuteForce();}  /// Go over all filaments in the network and call Calc. forces.
        
        for (auto it: (*_pmgdb) ) {it->CopmuteForce();} /// Go over all gost motors in the network and call Calc. forces.
        
        for (auto it: (*_pldb) ) {it->CopmuteForce();}  /// Go over all filaments in the network and call Calc. forces.
        
    }
    
    else{
        
        ResetForcesAux(*_pS->getBDB());
        for (auto it: (*_pfdb) ) {it->CopmuteForceAux();}
        
        for (auto it: (*_pmgdb) ) {it->CopmuteForceAux();}

        for (auto it: (*_pldb) ) {it->CopmuteForceAux();}
        
    }
    
    
}
