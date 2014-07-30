//
//  MSystem.cpp
//  CytoMech
//
//  Created by Konstantin Popov on 4/15/14.
//  Copyright (c) 2014 Konstantin Popov. All rights reserved.
//

#include "MSystem.h"
#include "MNetwork.h"

using namespace std;

System::System()
{
   	_pNetworkVector.push_back(AddNewNetwork());            //Create a Network (potentially can be more than one) upon initialization of a system;
	
}
//~System() {delete _pbdb, _pfdb;}



//void System::AddNewFilament(vector<double> v) { _pfilamentVector.push_back( _pfdb->CreateFilament(this, v));}


// Interfaces for creation of a Filament (vector with coordinates of beginig and end for the filament, and in what network (i) to create ):

void System::AddNewFilaments(vector<vector<vector<double> > > v, int i) {
    
    FilamentDB::Instance()->CreateFilament(this, _pNetworkVector[0], v);
    
}

/// Interfaces for creation of Linkers (vector with pointers to pairs of beads/one pair of beads, and in what network (i) to create ):
void System::AddNewLinkers(std::vector<std::vector<Cylinder *> > v, double streatchConst, int i){
    
    LinkerDB::Instance()->CreateLinker(this, _pNetworkVector[0], v, streatchConst);
    
}

//void System::AddNewLinker(Cylinder* pc1, Cylinder* pc2, double streatchConst, int i){
//
//    _pNetworkVector[i]->AddNewLinker(pc1, pc2, streatchConst);
//}

/// Interfaces to create a Ghost (no actual structuaral beads, only a force field between beads on filaments) (object that contains a vector with 4 beads, mechanical constants and methods for energy and force calculations):
void System::AddNewMotorGost(Cylinder* pc1, Cylinder* pc2, double k, double position1, double position2, int i){
    
    _pNetworkVector[i]->AddNewMotorGhost(pc1, pc2,k,position1, position2);

}

// For vector of motors(many motors):
void System::AddNewMotorGhosts(std::vector<std::vector<Cylinder* > > v, double k, double position1, double position2, int i){
    
    _pNetworkVector[i]->AddNewMotorGhosts(v,k,position1, position2);
}

/// Interface to remove a bead from a filament:
void System::RemoveBead(Bead*){
    
    
}

BeadDB* System::getBDB() {return _pbdb;}

CylinderDB* System::getCDB() {return _pcdb;}

int System::getSystemSize() {return _pbdb->size();}

Network* System::AddNewNetwork()
{
    Network* pn = new Network(this);
    return pn;
}

void System::CopmuteForce(int i)
{
    for (auto it: _pNetworkVector){
        
        (*it).CopmuteForce(i);
    }
}

double System::UpdateEnergy(double d)
{
    double u = 0;
    
    for (auto it: _pNetworkVector){
        
        u += (*it).CopmuteEnergy(d);
    }
    
    _energy = u;
    return u;
}

double System::getSystemEnergy() {return _energy;}
