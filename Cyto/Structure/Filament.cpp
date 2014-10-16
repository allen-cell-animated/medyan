//
//  Filament.cpp
//  CytoMech
//
//  Created by Konstantin Popov on 4/15/14.
//  Copyright (c) 2014 Konstantin Popov. All rights reserved.
//

#include "Filament.h"
#include "SystemParameters.h"
#include "Bead.h"
#include "BeadDB.h"
#include "CylinderDB.h"

using namespace std;
using namespace mathfunc;

Filament::Filament(SubSystem* ps, vector<double>& position, vector<double>& direction, int ID) : _ID(ID) {
   
    _pSubSystem = ps;
 
    ///Create beads
    Bead* b1 = BeadDB::Instance(BeadDBKey())->CreateBead(position, _beadIDPlusEnd++);
    auto pos2 = NextPointProjection(position, SystemParameters::Geometry().cylinderSize, direction);
    Bead* b2 = BeadDB::Instance(BeadDBKey())->CreateBead(pos2, _beadIDPlusEnd++);
    
    ///create cylinder
    Cylinder* c0 = CylinderDB::Instance(CylinderDBKey())->CreateCylinder(this, b1, b2);
    _pCylinderVector.push_back(c0);
    
    ///extend front
    auto newpoint = NextPointProjection(pos2, SystemParameters::Geometry().cylinderSize, direction);
    ExtendFront( newpoint );
    
#ifdef CHEMISTRY
    ///Update cylinder reactions
    for(auto &c : _pCylinderVector) { c->getCCylinder()->updateReactions(); }
#endif ///CHEMISTRY
}


Filament::Filament(SubSystem* ps, vector<vector<double> >& position, int numBeads, int ID, std::string projectionType) : _ID(ID) {

    _pSubSystem = ps;
    vector<vector<double> > tmpBeadsCoord;
    
    ///create a projection of beads
    if(projectionType == "STRAIGHT") tmpBeadsCoord = StraightFilamentProjection(position, numBeads);
    else if(projectionType == "ZIGZAG") tmpBeadsCoord = ZigZagFilamentProjection(position, numBeads);
    else {}
   
    ///Create beads
    
    auto direction = TwoPointDirection(tmpBeadsCoord[0], tmpBeadsCoord[1]);
    Bead* b1 = BeadDB::Instance(BeadDBKey())->CreateBead(tmpBeadsCoord[0], _beadIDPlusEnd++);
    Bead* b2 = BeadDB::Instance(BeadDBKey())->CreateBead(tmpBeadsCoord[1], _beadIDPlusEnd++);
    
    Cylinder* c0 = CylinderDB::Instance(CylinderDBKey())->CreateCylinder(this, b1, b2);
    _pCylinderVector.push_back(c0);
    
    for (int i = 2; i<numBeads; i++) {
        ExtendFront(tmpBeadsCoord[i]);  //Create n beads and n cylinders: x---x----x...x----x----o.
    }
    
#ifdef CHEMISTRY
    ///Update cylinder reactions
    for(auto &c : _pCylinderVector) { c->getCCylinder()->updateReactions(); }
#endif ///CHEMISTRY
}

///Extend front for initialization
void Filament::ExtendFront(vector<double>& coordinates) {
    
    Cylinder* cBack = _pCylinderVector.back();
    
    Bead* b2 = cBack->GetSecondBead();
    
    ///create a new bead
    auto direction = TwoPointDirection(b2->coordinate, coordinates);
    auto newBeadCoords = NextPointProjection(b2->coordinate, SystemParameters::Geometry().cylinderSize, direction);
    
    Bead* bNew = BeadDB::Instance(BeadDBKey())->CreateBead(newBeadCoords, _beadIDPlusEnd++);
    
    ///create cylinder
    Cylinder* c0 = CylinderDB::Instance(CylinderDBKey())->CreateCylinder(this, b2, bNew);
    c0->SetLast(true);
    _pCylinderVector.push_back(c0);
    
}

///Extend front for initialization
void Filament::ExtendBack(vector<double>& coordinates) {

    Cylinder* cFront = _pCylinderVector.front();
    Bead* b2 = cFront->GetFirstBead();
    
    //create a new bead
    auto direction = TwoPointDirection(b2->coordinate, coordinates);
    auto newBeadCoords = NextPointProjection(b2->coordinate, SystemParameters::Geometry().cylinderSize, direction);
    Bead* bNew = BeadDB::Instance(BeadDBKey())->CreateBead(newBeadCoords, _beadIDMinusEnd--);

    Cylinder* c0 = CylinderDB::Instance(CylinderDBKey())->CreateCylinder(this, bNew, b2);
    _pCylinderVector.push_front(c0);

}

///extend front at runtime
void Filament::ExtendFront() {
    if (_pCylinderVector.size()<1) {cout<<"ERROR FILAMENT TO SHORT"<<endl;}  /// Delete later
    
    else{
        Cylinder* cBack = _pCylinderVector.back();
        
        Bead* b1 = cBack->GetFirstBead();
        Bead* b2 = cBack->GetSecondBead();
        
        ///move last bead of last cylinder forward
        auto direction1 = TwoPointDirection(b1->coordinate, b2->coordinate);
        auto npp = NextPointProjection(b2->coordinate, SystemParameters::Geometry().monomerSize, direction1);
        
        ///create a new bead in same place as b2
        Bead* bNew = BeadDB::Instance(BeadDBKey())->CreateBead(npp, _beadIDPlusEnd++);
        
        Cylinder* c0 = CylinderDB::Instance(CylinderDBKey())->CreateCylinder(this, b2, bNew, true);
        _pCylinderVector.back()->SetLast(false);
        _pCylinderVector.push_back(c0);
        _pCylinderVector.back()->SetLast(true);
        
        _deltaPlusEnd++;
    }
}

///extend back at runtime
void Filament::ExtendBack() {
    if (_pCylinderVector.size()<1) {cout<<"ERROR FILAMENT TO SHORT"<<endl;}
    
    else{
        Cylinder* cFront = _pCylinderVector.front();
        
        Bead* b2 = cFront->GetFirstBead();
        Bead* b1 = cFront->GetSecondBead();
        
        ///move last bead of last cylinder forward
        auto direction1 = TwoPointDirection(b1->coordinate, b2->coordinate);
        auto npp = NextPointProjection(b2->coordinate, SystemParameters::Geometry().monomerSize, direction1);
        
        ///create a new bead in same place as b2
        Bead* bNew = BeadDB::Instance(BeadDBKey())->CreateBead(npp, _beadIDMinusEnd--);

        Cylinder* c0 = CylinderDB::Instance(CylinderDBKey())->CreateCylinder(this, bNew, b2,false, true);
        _pCylinderVector.push_front(c0);
        
        _deltaMinusEnd++;
    }
}

///Depolymerize front at runtime
void Filament::RetractFront() {
    if (_pCylinderVector.size()<1) {cout<<"ERROR FILAMENT TO SHORT"<<endl;}  /// Delete later
    
    else {
        
        Cylinder* retractionCylinder = _pCylinderVector.back();
        _pCylinderVector.pop_back();
        
        BeadDB::Instance(BeadDBKey())->RemoveBead(retractionCylinder->GetSecondBead());
        CylinderDB::Instance(CylinderDBKey())->RemoveCylinder(retractionCylinder);
        
        _pCylinderVector.back()->SetLast(true);

        _beadIDPlusEnd--;
        _deltaPlusEnd--;
    }
    
}

void Filament::RetractBack() {
    if (_pCylinderVector.size()<1) {cout<<"ERROR FILAMENT TO SHORT"<<endl;}  /// Delete later
    
    else {
        Cylinder* retractionCylinder = _pCylinderVector.front();
        _pCylinderVector.pop_front();
        
        BeadDB::Instance(BeadDBKey())->RemoveBead(retractionCylinder->GetFirstBead());
        CylinderDB::Instance(CylinderDBKey())->RemoveCylinder(retractionCylinder);
        
        _beadIDMinusEnd++;
        _deltaMinusEnd--;
    }
}

void Filament::PolymerizeFront() {
    
    Cylinder* cBack = _pCylinderVector.back();
    
    Bead* b1 = cBack->GetFirstBead();
    Bead* b2 = cBack->GetSecondBead();
    
    auto direction = TwoPointDirection(b1->coordinate, b2->coordinate);
    b2->coordinate = NextPointProjection(b2->coordinate, SystemParameters::Geometry().monomerSize, direction);
    
    ///increase length, update
#ifdef MECHANICS
    cBack->getMCylinder()->SetEqLength(cBack->getMCylinder()->GetEqLength()
                                       + SystemParameters::Geometry().monomerSize);
#endif
}

void Filament::PolymerizeBack() {
    
    Cylinder* cFront = _pCylinderVector.front();
    
    Bead* b2 = cFront->GetFirstBead();
    Bead* b1 = cFront->GetSecondBead();

    auto direction = TwoPointDirection(b1->coordinate, b2->coordinate);
    b2->coordinate = NextPointProjection(b2->coordinate, SystemParameters::Geometry().monomerSize, direction);

    ///increase length
#ifdef MECHANICS
    cFront->getMCylinder()->SetEqLength(cFront->getMCylinder()->GetEqLength()
                                        + SystemParameters::Geometry().monomerSize);
#endif
}

void Filament::DepolymerizeFront() {
    
    Cylinder* cBack = _pCylinderVector.back();
    
    Bead* b1 = cBack->GetFirstBead();
    Bead* b2 = cBack->GetSecondBead();

    auto direction = TwoPointDirection(b2->coordinate, b1->coordinate);
    b2->coordinate = NextPointProjection(b2->coordinate, SystemParameters::Geometry().monomerSize, direction);
    
    ///increase length, update
#ifdef MECHANICS
    cBack->getMCylinder()->SetEqLength(cBack->getMCylinder()->GetEqLength()
                                       - SystemParameters::Geometry().monomerSize);
#endif
}

void Filament::DepolymerizeBack() {
    
    Cylinder* cFront = _pCylinderVector.front();
    
    Bead* b2 = cFront->GetFirstBead();
    Bead* b1 = cFront->GetSecondBead();
    
    auto direction = TwoPointDirection(b2->coordinate, b1->coordinate);
    b2->coordinate = NextPointProjection(b2->coordinate, SystemParameters::Geometry().monomerSize, direction);
    
    ///increase length
#ifdef MECHANICS
    cFront->getMCylinder()->SetEqLength(cFront->getMCylinder()->GetEqLength()
                                        - SystemParameters::Geometry().monomerSize);
#endif
}


vector<vector<double> > Filament::StraightFilamentProjection(vector<vector<double>>& v, int numBeads){
    
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

vector<vector<double> > Filament::ZigZagFilamentProjection(vector<vector<double>>& v, int numBeads){
    
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



void Filament::DeleteBead(Bead*){
    
    cout<<"not implemented"<<endl;
}
