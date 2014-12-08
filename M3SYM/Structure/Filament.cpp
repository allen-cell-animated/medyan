
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

#include "Filament.h"

#include "Bead.h"
#include "Cylinder.h"

#include "SystemParameters.h"
#include "MathFunctions.h"

using namespace mathfunc;

Filament::Filament(SubSystem* s, vector<double>& position, vector<double>& direction) {
   
    _subSystem = s;
    
    //add to filament db
    FilamentDB::instance()->addFilament(this);
    _ID = FilamentDB::instance()->getFilamentID();
 
    //create beads
    Bead* b1 = new Bead(position, 0);
    auto pos2 = nextPointProjection(position, SystemParameters::Geometry().cylinderSize, direction);
    Bead* b2 = new Bead(pos2, 1);
    
    //create cylinder
    Cylinder* c0 = new Cylinder(this, b1, b2, 0);
    _cylinderVector.push_back(c0);
}


Filament::Filament(SubSystem* s, vector<vector<double> >& position, int numBeads, string projectionType){

    _subSystem = s;
    
    //add to filament db
    FilamentDB::instance()->addFilament(this);
    _ID = FilamentDB::instance()->getFilamentID();
    
    vector<vector<double> > tmpBeadsCoord;
    
    //create a projection of beads
    if(projectionType == "STRAIGHT") tmpBeadsCoord = straightFilamentProjection(position, numBeads);
    else if(projectionType == "ZIGZAG") tmpBeadsCoord = zigZagFilamentProjection(position, numBeads);
    
    else {}
   
    //create beads
    auto direction = twoPointDirection(tmpBeadsCoord[0], tmpBeadsCoord[1]);
    Bead* b1 = new Bead(tmpBeadsCoord[0], 0);
    Bead* b2 = new Bead(tmpBeadsCoord[1], 1);
    
    Cylinder* c0 = new Cylinder(this, b1, b2, 0);
    _cylinderVector.push_back(c0);
    
    for (int i = 2; i<numBeads; i++) {
        extendFront(tmpBeadsCoord[i]);  //Create n beads and n cylinders: x---x----x...x----x----o.
    }
}

Filament::~Filament() {
    
    //remove cylinders, beads from system
    for(auto &c : _cylinderVector) {
        //remove first bead
        delete c->getFirstBead();
        //remove second bead if last
        if(c->last()) delete c->getSecondBead();
        //remove cylinder
        delete c;
    }
}


//Extend front for initialization
void Filament::extendFront(vector<double>& coordinates) {
    
    Cylinder* cBack = _cylinderVector.back();
    Bead* b2 = cBack->getSecondBead();
    
    //create a new bead
    auto direction = twoPointDirection(b2->coordinate, coordinates);
    auto newBeadCoords = nextPointProjection(b2->coordinate, SystemParameters::Geometry().cylinderSize, direction);
    
    Bead* bNew = new Bead(newBeadCoords, b2->getPositionFilament() + 1);
    
    //create cylinder
    Cylinder* c0 = new Cylinder(this, b2, bNew, _cylinderVector.size());
    c0->setLast(true);
    _cylinderVector.push_back(c0);
    
}

//Extend front for initialization
void Filament::extendBack(vector<double>& coordinates) {

    Cylinder* cFront = _cylinderVector.front();
    int lastPositionFilament = cFront->getPositionFilament();
    Bead* b2 = cFront->getFirstBead();
    
    //create a new bead
    auto direction = twoPointDirection(b2->coordinate, coordinates);
    auto newBeadCoords = nextPointProjection(b2->coordinate, SystemParameters::Geometry().cylinderSize, direction);
    Bead* bNew = new Bead(newBeadCoords, b2->getPositionFilament() - 1);
    
    Cylinder* c0 = new Cylinder(this, bNew, b2, lastPositionFilament - 1);
    _cylinderVector.push_front(c0);

}

//extend front at runtime
void Filament::extendFront() {
    if (_cylinderVector.size()<1) {cout<<"ERROR FILAMENT TO SHORT"<<endl;}  // Delete later
    
    else{
        Cylinder* cBack = _cylinderVector.back();
        
        Bead* b1 = cBack->getFirstBead();
        Bead* b2 = cBack->getSecondBead();
        
        //move last bead of last cylinder forward
        auto direction1 = twoPointDirection(b1->coordinate, b2->coordinate);
        auto npp = nextPointProjection(b2->coordinate, SystemParameters::Geometry().monomerSize, direction1);
        
        //create a new bead in same place as b2
        Bead* bNew = new Bead(npp, b2->getPositionFilament() + 1);
        
        Cylinder* c0 = new Cylinder(this, b2, bNew, _cylinderVector.size(), true);
        _cylinderVector.back()->setLast(false);
        _cylinderVector.push_back(c0);
        _cylinderVector.back()->setLast(true);
        
        _deltaPlusEnd++;
    }
}

//extend back at runtime
void Filament::extendBack() {
    if (_cylinderVector.size()<1) {cout<<"ERROR FILAMENT TO SHORT"<<endl;}
    
    else{
        Cylinder* cFront = _cylinderVector.front();
        int lastPositionFilament = cFront->getPositionFilament();
        
        Bead* b2 = cFront->getFirstBead();
        Bead* b1 = cFront->getSecondBead();
        
        //move last bead of last cylinder forward
        auto direction1 = twoPointDirection(b1->coordinate, b2->coordinate);
        auto npp = nextPointProjection(b2->coordinate, SystemParameters::Geometry().monomerSize, direction1);
        
        //create a new bead in same place as b2
        Bead* bNew = new Bead(npp, b2->getPositionFilament() - 1);

        Cylinder* c0 = new Cylinder(this, bNew, b2, lastPositionFilament - 1, false, true);
        _cylinderVector.push_front(c0);
        
        _deltaMinusEnd++;
    }
}

//Depolymerize front at runtime
void Filament::retractFront() {
    if (_cylinderVector.size()<1) {cout<<"ERROR FILAMENT TO SHORT"<<endl;}  // Delete later
    
    else {
        Cylinder* retractionCylinder = _cylinderVector.back();
        _cylinderVector.pop_back();
        
        delete retractionCylinder->getSecondBead();
        delete retractionCylinder;
        
        _cylinderVector.back()->setLast(true);

        _deltaPlusEnd--;
    }
    
}

void Filament::retractBack() {
    if (_cylinderVector.size()<1) {cout<<"ERROR FILAMENT TO SHORT"<<endl;}  // Delete later
    
    else {
        Cylinder* retractionCylinder = _cylinderVector.front();
        _cylinderVector.pop_front();
        
        delete retractionCylinder->getFirstBead();
        delete retractionCylinder;
        
        _deltaMinusEnd--;
    }
}

void Filament::polymerizeFront() {
    
    Cylinder* cBack = _cylinderVector.back();
    
    Bead* b1 = cBack->getFirstBead();
    Bead* b2 = cBack->getSecondBead();
    
    auto direction = twoPointDirection(b1->coordinate, b2->coordinate);
    b2->coordinate = nextPointProjection(b2->coordinate, SystemParameters::Geometry().monomerSize, direction);
    
    //increase length, update
#ifdef MECHANICS
    cBack->getMCylinder()->setEqLength(cBack->getMCylinder()->getEqLength() + SystemParameters::Geometry().monomerSize);
#endif
}

void Filament::polymerizeBack() {
    
    Cylinder* cFront = _cylinderVector.front();
    
    Bead* b2 = cFront->getFirstBead();
    Bead* b1 = cFront->getSecondBead();

    auto direction = twoPointDirection(b1->coordinate, b2->coordinate);
    b2->coordinate = nextPointProjection(b2->coordinate, SystemParameters::Geometry().monomerSize, direction);

    //increase length
#ifdef MECHANICS
    cFront->getMCylinder()->setEqLength(cFront->getMCylinder()->getEqLength() + SystemParameters::Geometry().monomerSize);
#endif
}

void Filament::depolymerizeFront() {
    
    Cylinder* cBack = _cylinderVector.back();
    
    Bead* b1 = cBack->getFirstBead();
    Bead* b2 = cBack->getSecondBead();

    auto direction = twoPointDirection(b2->coordinate, b1->coordinate);
    b2->coordinate = nextPointProjection(b2->coordinate, SystemParameters::Geometry().monomerSize, direction);
    
    //increase length, update
#ifdef MECHANICS
    cBack->getMCylinder()->setEqLength(cBack->getMCylinder()->getEqLength() - SystemParameters::Geometry().monomerSize);
#endif
}

void Filament::depolymerizeBack() {
    
    Cylinder* cFront = _cylinderVector.front();
    
    Bead* b2 = cFront->getFirstBead();
    Bead* b1 = cFront->getSecondBead();
    
    auto direction = twoPointDirection(b2->coordinate, b1->coordinate);
    b2->coordinate = nextPointProjection(b2->coordinate, SystemParameters::Geometry().monomerSize, direction);
    
    //increase length
#ifdef MECHANICS
    cFront->getMCylinder()->setEqLength(cFront->getMCylinder()->getEqLength() - SystemParameters::Geometry().monomerSize);
#endif
}

void Filament::printChemComposition() {
    for (auto &c : _cylinderVector) {
        c->getCCylinder()->printCCylinder();
    }
}

vector<vector<double> > Filament::straightFilamentProjection(vector<vector<double>>& v, int numBeads){
    
    vector<vector<double>> coordinate;
    vector<double> tmpVec (3, 0);
    vector<double> tau (3, 0);
    double invD = 1/twoPointDistance(v[1], v[0]);
    tau[0] = invD * ( v[1][0] - v[0][0] );
    tau[1] = invD * ( v[1][1] - v[0][1] );
    tau[2] = invD * ( v[1][2] - v[0][2] );
    
    for (int i = 0; i<numBeads; i++) {
        
        tmpVec[0] = v[0][0] + SystemParameters::Geometry().cylinderSize * i * tau[0];
        tmpVec[1] = v[0][1] + SystemParameters::Geometry().cylinderSize * i * tau[1];
        tmpVec[2] = v[0][2] + SystemParameters::Geometry().cylinderSize * i * tau[2];
        
        coordinate.push_back(tmpVec);
    }
    return coordinate;
}

vector<vector<double> > Filament::zigZagFilamentProjection(vector<vector<double>>& v, int numBeads){
    
    vector<vector<double>> coordinate;
    vector<double> tmpVec (3, 0);
    vector<double> tau (3, 0);
    double invD = 1/twoPointDistance(v[1], v[0]);
    tau[0] = invD * ( v[1][0] - v[0][0] );
    tau[1] = invD * ( v[1][1] - v[0][1] );
    tau[2] = invD * ( v[1][2] - v[0][2] );
    
    vector<double> perptau = {-tau[1], tau[0], tau[2]};
    
    
    for (int i = 0; i<numBeads; i++) {
        
        if(i%2 == 0) {
            tmpVec[0] = v[0][0] + SystemParameters::Geometry().cylinderSize * i * tau[0];
            tmpVec[1] = v[0][1] + SystemParameters::Geometry().cylinderSize * i * tau[1];
            tmpVec[2] = v[0][2] + SystemParameters::Geometry().cylinderSize * i * tau[2];
        }
        else {
            tmpVec[0] = v[0][0] + SystemParameters::Geometry().cylinderSize * i * perptau[0];
            tmpVec[1] = v[0][1] + SystemParameters::Geometry().cylinderSize * i * perptau[1];
            tmpVec[2] = v[0][2] + SystemParameters::Geometry().cylinderSize * i * perptau[2];
        }
        
        coordinate.push_back(tmpVec);
    }
    return coordinate;
}

