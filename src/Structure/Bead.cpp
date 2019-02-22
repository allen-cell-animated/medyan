
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v3.1
//
//  Copyright (2015-2016)  Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the MEDYAN web page for more information:
//  http://www.medyan.org
//------------------------------------------------------------------

#include "Bead.h"

#include "Compartment.h"
#include "Composite.h"

#include "SysParams.h"
#include "GController.h"
#include "MathFunctions.h"

using namespace mathfunc;

Database<Bead*> Bead::_beads;
Database<Bead*> Bead::_pinnedBeads;
//static vars needed to vectorize on-the-fly
int Bead::maxbindex = 0;
int Bead::vectormaxsize = 0;
int Bead::Nbeads = 0;
bool Bead::triggercylindervectorization = false;
vector<int> Bead::removedbindex;//vector of bead indices that were once alloted to other
// beads but are free to be reallocated now.

Bead::Bead (vector<double> v, Composite* parent, int position)
// Qin add brforce, pinforce
    : Trackable(true),
      coordinate(v), coordinateP(v),
      force(3, 0), forceAux(3, 0), forceAuxP(3, 0), brforce(3, 0), pinforce(3,0),
      _position(position), _birthTime(tau()),_ID(_beads.getID()) {
    
    parent->addChild(unique_ptr<Component>(this));
          
    loadForcesP = vector<double>(SysParams::Geometry().cylinderNumMon[getType()], 0.0);
    loadForcesM = vector<double>(SysParams::Geometry().cylinderNumMon[getType()], 0.0);
    
    //Find compartment
    try {_compartment = GController::getCompartment(v);}
    catch (exception& e) {
        
        cout << e.what() << endl;
        
        printSelf();
        
        //also print parent info
        getParent()->printSelf();
        
        exit(EXIT_FAILURE);
    }

    //revectorize if needed
    revectorizeifneeded();
    //set bindex
    if(removedbindex.size() == 0)
    {_dbIndex = maxbindex;
    maxbindex++;
    }
    else{
//        std::cout<<"reusing bIndex "<<removedbindex.at(0)<<endl;
        _dbIndex = removedbindex.at(0);
        removedbindex.erase(removedbindex.begin());
    }
    Nbeads = _beads.getElements().size();

    //copy bead coordiantes to the appropriate spot in the coord vector.
    copycoordinatestovector();
}

Bead::Bead(Composite* parent, int position)
// Qin add brforce, pinforce
    : Trackable(true),
    coordinate(3, 0), coordinateP(3, 0),
    force(3, 0), forceAux(3, 0), forceAuxP(3, 0), brforce(3, 0), pinforce(3,0), _position(position) {
    
    parent->addChild(unique_ptr<Component>(this));
    //check if you need to revectorize.
    revectorizeifneeded();
    //set bindex based on maxbindex if there were no beads removed.
    if(removedbindex.size() == 0)
    {_dbIndex = maxbindex;
        maxbindex++;
    }
    // if beads were removed earlier, allot one of the available bead indices.
    else{
//        std::cout<<"reusing bIndex v2"<<removedbindex.at(0)<<endl;
        _dbIndex = removedbindex.at(0);
        removedbindex.erase(removedbindex.begin());
    }
    Nbeads = _beads.getElements().size();
    //copy bead coordiantes to the appropriate spot in the coord vector.
    copycoordinatestovector();
}

void Bead::updatePosition() {
    
    try {GController::getCompartment(coordinate);}
    catch (exception& e) {
        
        //print exception
        cout << e.what() << endl;
        
        printSelf();
        
        //also print parent info
        getParent()->printSelf();
        
        //exit
        exit(EXIT_FAILURE);
    }
}

void Bead::printSelf() {
    
    cout << endl;
    
    cout << "Bead: ptr = " << this << endl;
    cout << "Coordinates = " << coordinate[0] << ", " << coordinate[1] << ", " << coordinate[2] << endl;
    cout << "Previous coordinates before minimization = " << coordinateP[0] << ", " << coordinateP[1] << ", " << coordinateP[2] << endl;
    cout << "Forces = " << force[0] << ", " << force[1] << ", " << force[2] << endl;
    cout << "Auxiliary forces = " << forceAux[0] << ", " << forceAux[1] << ", " << forceAux[2] << endl;

    cout << "Position on structure = " << _position << endl;
    cout << "Birth time = " << _birthTime << endl;
    
    
    cout << endl;
}

double Bead::getLoadForcesP() {
    
    if (lfip < 0)
        return loadForcesP[0];
        
    if (lfip >= loadForcesP.size())
        return loadForcesP.back();
    
    else return loadForcesP[lfip];
}

double Bead::getLoadForcesM() {
    
    if (lfim < 0)
        return loadForcesM[0];
    
    if (lfim >= loadForcesM.size())
        return loadForcesM.back();
    
    else return loadForcesM[lfim];
}

