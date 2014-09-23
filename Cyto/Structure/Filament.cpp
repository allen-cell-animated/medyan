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

Filament::Filament(SubSystem* ps, vector<double>& position, vector<double>& direction){
   
    _pSubSystem = ps;
 
    ///Create cylinders, beads
    Bead* b1 = BeadDB::Instance(BeadDBKey())->CreateBead(position);
    auto npp = NextPointProjection(b1->coordinate, SystemParameters::Geometry().cylinderSize, direction);
    auto midpoint = MidPointCoordinate(b1->coordinate, npp, 0.5);
    Compartment* c;
    try {c =  GController::getCompartment(midpoint);}
    catch (exception& e) {std:: cout << e.what(); exit(EXIT_FAILURE);}
    
    Cylinder* c1 = CylinderDB::Instance(CylinderDBKey())->CreateCylinder(this, b1, c);
    _pLastCylinder = c1;
    auto point = NextPointProjection(position, SystemParameters::Geometry().cylinderSize, direction);
    PolymerizeFront( point );
}


Filament::Filament(SubSystem* ps, vector<vector<double> >& position, int numBeads){

    _pSubSystem = ps;
    
    vector<vector<double> > tmpBeadsCoord = StraightFilamentProjection(position, numBeads);
    //this function calculate coordinates for all beads on the line separated by a segment length.
   
    ///Create beads, cylinders
    Bead* b0 = BeadDB::Instance(BeadDBKey())->CreateBead(tmpBeadsCoord[0]);
    
    auto midpoint = MidPointCoordinate(tmpBeadsCoord[0], tmpBeadsCoord[1], 0.5);
    Compartment* c;
    try {c = GController::getCompartment(midpoint);}
    catch (exception& e) {std:: cout << e.what(); exit(EXIT_FAILURE);}
    
    Cylinder* c0 = CylinderDB::Instance(CylinderDBKey())->CreateCylinder(this, b0, c);
    _pLastCylinder = c0;
   
    
    for (int i = 1; i<numBeads; i++) {
        PolymerizeFront(tmpBeadsCoord[i]);  //Create n beads and n cylinders: x---x----x...x----x----o.
    }
    
    
}

///Polymerize front for initialization
void Filament::PolymerizeFront(vector<double>& coordinates) {
    
    Bead* b = BeadDB::Instance(BeadDBKey())->CreateBead(coordinates);

    auto tau = TwoPointDirection( _pLastCylinder->getMCylinder()->GetFirstBead()->coordinate, coordinates);
    auto npp = NextPointProjection(coordinates, SystemParameters::Geometry().cylinderSize, tau);
    auto midpoint = MidPointCoordinate(coordinates, npp, 0.5);

    Compartment* c;
    try {c = GController::getCompartment(midpoint);}
    catch (exception& e) {std:: cout << e.what(); exit(EXIT_FAILURE);}
    
    Cylinder* c0 = CylinderDB::Instance(CylinderDBKey())->CreateCylinder(this, b, c);
    _pLastCylinder->getMCylinder()->SetSecondBead(b);
    _pLastCylinder->SetLast(false);
    
    _pCylinderVector.push_back(_pLastCylinder);
    
    _pLastCylinder = c0;
    
}

///Polymerize front for initialization
void Filament::PolymerizeBack(vector<double>& coordinates) {
    
    Bead* b = BeadDB::Instance(BeadDBKey())->CreateBead(coordinates);
    
    auto midpoint = MidPointCoordinate(coordinates, _pCylinderVector[0]->getMCylinder()->GetFirstBead()->coordinate, 0.5);
    Compartment* c;
    try {c = GController::getCompartment(midpoint);}
    catch (exception& e) {std:: cout << e.what(); exit(EXIT_FAILURE);}
    
    Cylinder* c0 = CylinderDB::Instance(CylinderDBKey())->CreateCylinder(this, b, c);
    c0->getMCylinder()->SetSecondBead(_pCylinderVector[0]->getMCylinder()->GetFirstBead());
    c0->SetLast(false);
    _pCylinderVector.push_front(c0);

}

///Polymerize front at runtime
void Filament::PolymerizeFront() {
    if (_pCylinderVector.size()<1) {cout<<"ERROR FILAMENT TO SHORT"<<endl;}  /// Delete later
    
    else{
        auto tau = TwoPointDirection( _pCylinderVector[_pCylinderVector.size() - 1]->getMCylinder()->GetFirstBead()->coordinate,
                                      _pCylinderVector[_pCylinderVector.size() - 1]->getMCylinder()->GetSecondBead()->coordinate);
        auto npp1 = NextPointProjection(
                    _pCylinderVector[_pCylinderVector.size()-1]->getMCylinder()->GetSecondBead()->coordinate,
                    SystemParameters::Geometry().cylinderSize, tau);
        
        Bead* b = BeadDB::Instance(BeadDBKey())->CreateBead(npp1);
        
        auto npp2 = NextPointProjection(b->coordinate, SystemParameters::Geometry().cylinderSize, tau);
        auto midpoint = MidPointCoordinate(b->coordinate, npp2, 0.5);
        Compartment* c;
        
        try {c = GController::getCompartment(midpoint);}
        catch (exception& e) {std:: cout << e.what(); exit(EXIT_FAILURE);}
        
        Cylinder* c0 = CylinderDB::Instance(CylinderDBKey())->CreateCylinder(this, b, c, true);
        _pLastCylinder->getMCylinder()->SetSecondBead(b);
        _pLastCylinder->SetLast(false);
        _pCylinderVector.push_back(_pLastCylinder);
        _pLastCylinder = c0;
                                       
    }
}

///Polymerize back at runtime
void Filament::PolymerizeBack() {
    if (_pCylinderVector.size()<1) {cout<<"ERROR FILAMENT TO SHORT"<<endl;}
    
    else{
        /// Check out!!!!
        auto tau = TwoPointDirection(
         _pCylinderVector[0]->getMCylinder()->GetSecondBead()->coordinate,
         _pCylinderVector[0]->getMCylinder()->GetFirstBead()->coordinate);

        auto npp = NextPointProjection(_pCylinderVector[0]->getMCylinder()->GetFirstBead()->coordinate,
                                       SystemParameters::Geometry().cylinderSize, tau);
        Bead* b = BeadDB::Instance(BeadDBKey())->CreateBead(npp);
        
        auto midpoint = MidPointCoordinate(b->coordinate, _pCylinderVector[0]->getMCylinder()->GetFirstBead()->coordinate, 0.5);
        Compartment* c;
        try {c = GController::getCompartment(midpoint);}
        catch (exception& e) {std:: cout << e.what(); exit(EXIT_FAILURE);}
        
        Cylinder* c0 = CylinderDB::Instance(CylinderDBKey())->CreateCylinder(this, b, c, false, true);
        c0->getMCylinder()->SetSecondBead(_pCylinderVector[0]->getMCylinder()->GetFirstBead());
        c0->SetLast(false);
        _pCylinderVector.push_front(c0);
    }
}

///Depolymerize front at runtime
void Filament::DepolymerizeFront() {
    if (_pCylinderVector.size()<1) {cout<<"ERROR FILAMENT TO SHORT"<<endl;}  /// Delete later
    
    else {
        CylinderDB::Instance(CylinderDBKey())->RemoveCylinder(_pLastCylinder);
        BeadDB::Instance(BeadDBKey())->RemoveBead(_pCylinderVector[_pCylinderVector.size() - 1]->getMCylinder()->GetSecondBead());
        
        _pLastCylinder = _pCylinderVector[_pCylinderVector.size() - 1];
        _pLastCylinder->getMCylinder()->SetSecondBead(nullptr);
        _pLastCylinder->SetLast(true);
        _pCylinderVector.pop_back();
    }
    
}

void Filament::DepolymerizeBack() {
    if (_pCylinderVector.size()<1) {cout<<"ERROR FILAMENT TO SHORT"<<endl;}  /// Delete later
    
    else {
        BeadDB::Instance(BeadDBKey())->RemoveBead(_pCylinderVector[0]->getMCylinder()->GetSecondBead());
        CylinderDB::Instance(CylinderDBKey())->RemoveCylinder(_pCylinderVector[0]);
        
        _pCylinderVector.pop_front();
    }
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


void Filament::DeleteBead(Bead*){
    
    cout<<"not implemented"<<endl;
}
