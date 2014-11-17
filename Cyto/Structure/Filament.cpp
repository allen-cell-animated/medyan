//
//  Filament.cpp
//  CytoMech
//
//  Created by Konstantin Popov on 4/15/14.
//  Copyright (c) 2014 Konstantin Popov. All rights reserved.
//

#include "Filament.h"

#include "Bead.h"
#include "BeadDB.h"
#include "CylinderDB.h"

#include "SystemParameters.h"
#include "MathFunctions.h"

using namespace mathfunc;

Filament::Filament(SubSystem* s, vector<double>& position, vector<double>& direction, int ID) : _ID(ID) {
   
    _subSystem = s;
 
    ///Create beads
    Bead* b1 = BeadDB::instance(BeadDBKey())->createBead(position, 0);
    auto pos2 = NextPointProjection(position, SystemParameters::Geometry().cylinderSize, direction);
    Bead* b2 = BeadDB::instance(BeadDBKey())->createBead(pos2, 1);
    
    ///create cylinder
    Cylinder* c0 = CylinderDB::instance(CylinderDBKey())->createCylinder(this, b1, b2, 0);
    _cylinderVector.push_back(c0);
}


Filament::Filament(SubSystem* s, vector<vector<double> >& position, int numBeads, int ID, string projectionType) : _ID(ID) {

    _subSystem = s;
    vector<vector<double> > tmpBeadsCoord;
    
    ///create a projection of beads
    if(projectionType == "STRAIGHT") tmpBeadsCoord = straightFilamentProjection(position, numBeads);
    else if(projectionType == "ZIGZAG") tmpBeadsCoord = zigZagFilamentProjection(position, numBeads);
    else {}
   
    ///Create beads
    auto direction = TwoPointDirection(tmpBeadsCoord[0], tmpBeadsCoord[1]);
    Bead* b1 = BeadDB::instance(BeadDBKey())->createBead(tmpBeadsCoord[0], 0);
    Bead* b2 = BeadDB::instance(BeadDBKey())->createBead(tmpBeadsCoord[1], 1);
    
    Cylinder* c0 = CylinderDB::instance(CylinderDBKey())->createCylinder(this, b1, b2, 0);
    _cylinderVector.push_back(c0);
    
    for (int i = 2; i<numBeads; i++) {
        extendFront(tmpBeadsCoord[i]);  //Create n beads and n cylinders: x---x----x...x----x----o.
    }
}

Filament::~Filament() {
    
    ///remove cylinders, beads from system
    for(auto &c : _cylinderVector) {
        //remove first bead
        BeadDB::instance(BeadDBKey())->removeBead(c->getFirstBead());
        //remove second bead if last
        if(c->last()) BeadDB::instance(BeadDBKey())->removeBead(c->getSecondBead());
        //remove cylinder
        CylinderDB::instance(CylinderDBKey())->removeCylinder(c);
    }
}


///Extend front for initialization
void Filament::extendFront(vector<double>& coordinates) {
    
    Cylinder* cBack = _cylinderVector.back();
    Bead* b2 = cBack->getSecondBead();
    
    ///create a new bead
    auto direction = TwoPointDirection(b2->coordinate, coordinates);
    auto newBeadCoords = NextPointProjection(b2->coordinate, SystemParameters::Geometry().cylinderSize, direction);
    
    Bead* bNew = BeadDB::instance(BeadDBKey())->createBead(newBeadCoords, b2->getPositionFilament() + 1);
    
    ///create cylinder
    Cylinder* c0 = CylinderDB::instance(CylinderDBKey())->createCylinder(this, b2, bNew, _cylinderVector.size());
    c0->setLast(true);
    _cylinderVector.push_back(c0);
    
}

///Extend front for initialization
void Filament::extendBack(vector<double>& coordinates) {

    Cylinder* cFront = _cylinderVector.front();
    int lastPositionFilament = cFront->getPositionFilament();
    Bead* b2 = cFront->getFirstBead();
    
    //create a new bead
    auto direction = TwoPointDirection(b2->coordinate, coordinates);
    auto newBeadCoords = NextPointProjection(b2->coordinate, SystemParameters::Geometry().cylinderSize, direction);
    Bead* bNew = BeadDB::instance(BeadDBKey())->createBead(newBeadCoords, b2->getPositionFilament() - 1);
    
    Cylinder* c0 = CylinderDB::instance(CylinderDBKey())->createCylinder(this, bNew, b2, lastPositionFilament - 1);
    _cylinderVector.push_front(c0);

}

///extend front at runtime
void Filament::extendFront() {
    if (_cylinderVector.size()<1) {cout<<"ERROR FILAMENT TO SHORT"<<endl;}  /// Delete later
    
    else{
        Cylinder* cBack = _cylinderVector.back();
        
        Bead* b1 = cBack->getFirstBead();
        Bead* b2 = cBack->getSecondBead();
        
        ///move last bead of last cylinder forward
        auto direction1 = TwoPointDirection(b1->coordinate, b2->coordinate);
        auto npp = NextPointProjection(b2->coordinate, SystemParameters::Geometry().monomerSize, direction1);
        
        ///create a new bead in same place as b2
        Bead* bNew = BeadDB::instance(BeadDBKey())->createBead(npp, b2->getPositionFilament() + 1);
        
        Cylinder* c0 = CylinderDB::instance(CylinderDBKey())->createCylinder(this, b2, bNew, _cylinderVector.size(), true);
        _cylinderVector.back()->setLast(false);
        _cylinderVector.push_back(c0);
        _cylinderVector.back()->setLast(true);
        
        _deltaPlusEnd++;
    }
}

///extend back at runtime
void Filament::extendBack() {
    if (_cylinderVector.size()<1) {cout<<"ERROR FILAMENT TO SHORT"<<endl;}
    
    else{
        Cylinder* cFront = _cylinderVector.front();
        int lastPositionFilament = cFront->getPositionFilament();
        
        Bead* b2 = cFront->getFirstBead();
        Bead* b1 = cFront->getSecondBead();
        
        ///move last bead of last cylinder forward
        auto direction1 = TwoPointDirection(b1->coordinate, b2->coordinate);
        auto npp = NextPointProjection(b2->coordinate, SystemParameters::Geometry().monomerSize, direction1);
        
        ///create a new bead in same place as b2
        Bead* bNew = BeadDB::instance(BeadDBKey())->createBead(npp, b2->getPositionFilament() - 1);

        Cylinder* c0 = CylinderDB::instance(CylinderDBKey())->createCylinder(this, bNew, b2, lastPositionFilament - 1, false, true);
        _cylinderVector.push_front(c0);
        
        _deltaMinusEnd++;
    }
}

///Depolymerize front at runtime
void Filament::retractFront() {
    if (_cylinderVector.size()<1) {cout<<"ERROR FILAMENT TO SHORT"<<endl;}  /// Delete later
    
    else {
        
        Cylinder* retractionCylinder = _cylinderVector.back();
        _cylinderVector.pop_back();
        
        BeadDB::instance(BeadDBKey())->removeBead(retractionCylinder->getSecondBead());
        CylinderDB::instance(CylinderDBKey())->removeCylinder(retractionCylinder);
        
        _cylinderVector.back()->setLast(true);

        _deltaPlusEnd--;
    }
    
}

void Filament::retractBack() {
    if (_cylinderVector.size()<1) {cout<<"ERROR FILAMENT TO SHORT"<<endl;}  /// Delete later
    
    else {
        Cylinder* retractionCylinder = _cylinderVector.front();
        _cylinderVector.pop_front();
        
        BeadDB::instance(BeadDBKey())->removeBead(retractionCylinder->getFirstBead());
        CylinderDB::instance(CylinderDBKey())->removeCylinder(retractionCylinder);
        
        _deltaMinusEnd--;
    }
}

void Filament::polymerizeFront() {
    
    Cylinder* cBack = _cylinderVector.back();
    
    Bead* b1 = cBack->getFirstBead();
    Bead* b2 = cBack->getSecondBead();
    
    auto direction = TwoPointDirection(b1->coordinate, b2->coordinate);
    b2->coordinate = NextPointProjection(b2->coordinate, SystemParameters::Geometry().monomerSize, direction);
    
    ///increase length, update
#ifdef MECHANICS
    cBack->getMCylinder()->setEqLength(cBack->getMCylinder()->getEqLength() + SystemParameters::Geometry().monomerSize);
#endif
}

void Filament::polymerizeBack() {
    
    Cylinder* cFront = _cylinderVector.front();
    
    Bead* b2 = cFront->getFirstBead();
    Bead* b1 = cFront->getSecondBead();

    auto direction = TwoPointDirection(b1->coordinate, b2->coordinate);
    b2->coordinate = NextPointProjection(b2->coordinate, SystemParameters::Geometry().monomerSize, direction);

    ///increase length
#ifdef MECHANICS
    cFront->getMCylinder()->setEqLength(cFront->getMCylinder()->getEqLength() + SystemParameters::Geometry().monomerSize);
#endif
}

void Filament::depolymerizeFront() {
    
    Cylinder* cBack = _cylinderVector.back();
    
    Bead* b1 = cBack->getFirstBead();
    Bead* b2 = cBack->getSecondBead();

    auto direction = TwoPointDirection(b2->coordinate, b1->coordinate);
    b2->coordinate = NextPointProjection(b2->coordinate, SystemParameters::Geometry().monomerSize, direction);
    
    ///increase length, update
#ifdef MECHANICS
    cBack->getMCylinder()->setEqLength(cBack->getMCylinder()->getEqLength() - SystemParameters::Geometry().monomerSize);
#endif
}

void Filament::depolymerizeBack() {
    
    Cylinder* cFront = _cylinderVector.front();
    
    Bead* b2 = cFront->getFirstBead();
    Bead* b1 = cFront->getSecondBead();
    
    auto direction = TwoPointDirection(b2->coordinate, b1->coordinate);
    b2->coordinate = NextPointProjection(b2->coordinate, SystemParameters::Geometry().monomerSize, direction);
    
    ///increase length
#ifdef MECHANICS
    cFront->getMCylinder()->setEqLength(cFront->getMCylinder()->getEqLength() - SystemParameters::Geometry().monomerSize);
#endif
}


vector<vector<double> > Filament::straightFilamentProjection(vector<vector<double>>& v, int numBeads){
    
    vector<vector<double>> coordinate;
    vector<double> tmpVec (3, 0);
    vector<double> tau (3, 0);
    double invD = 1/TwoPointDistance(v[1], v[0]);
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
    double invD = 1/TwoPointDistance(v[1], v[0]);
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



void Filament::deleteBead(Bead*){
    
    cout<<"not implemented"<<endl;
}
