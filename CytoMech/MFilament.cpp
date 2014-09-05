//
//  MFilament.cpp
//  CytoMech
//
//  Created by Konstantin Popov on 4/15/14.
//  Copyright (c) 2014 Konstantin Popov. All rights reserved.
//

#include "MFilament.h"
#include "MBead.h"
#include "SubSystem.h"

using namespace std;
using namespace mathfunc;



Filament::Filament(SubSystem* ps, vector<double> position, vector<double> direction){
   
    
    //This constructor creates a short filament, containing only two beads. Coordinates of the first bead is an input, second is set up by using an input direction and a coarsegrain length L. Using all this, two constructors for beads are called.
    
    _pSubSystem = ps;
 
	
    
    Bead* b1 = BeadDB::Instance(BeadDBKey())->CreateBead(position);
    Cylinder* c1 = CylinderDB::Instance(CylinderDBKey())->CreateCylinder(this, b1);
    _pLastCylinder = c1;
   

    
    PolymerizeFront( NextPointProjection(position, L, direction) );
}


Filament::Filament(SubSystem* ps, vector<vector<double> > position, int numBeads){
    
    /// This constructor is called to create a longer filament. It creates a filament with a number of beads numBeads. Filaments starts and ends in the point determined by position vector and has a direction direction. Number of beads is equal to the number of cylinders. The last cylinder doesnt have an end(second) bead and will not be pushed to cylinder vector, but will be stored in the _pLastCylinder;
    
    _pSubSystem = ps;
   
    
    vector<vector<double> > tmpBeadsCoord = StraightFilamentProjection(position, numBeads); //this function calculate coordinates for all beads on the line separated by a segment length.
   
    Bead* b0 = BeadDB::Instance(BeadDBKey())->CreateBead(tmpBeadsCoord[0]);
    Cylinder* c0 = CylinderDB::Instance(CylinderDBKey())->CreateCylinder(this, b0);
    _pLastCylinder = c0;
   
    
    for (int i = 1; i<numBeads; i++) {
        
        PolymerizeFront(tmpBeadsCoord[i]);  //Create n beads and n cylinders: x---x----x...x----x----o.
        }
    
}

void Filament::PolymerizeFront(vector<double> coordinates) {
    
    Bead* b = BeadDB::Instance(BeadDBKey())->CreateBead(coordinates);
    Cylinder* c = CylinderDB::Instance(CylinderDBKey())->CreateCylinder(this, b);
    _pLastCylinder->getMCylinder()->SetSecondBead(b);
    _pLastCylinder->SetLast(false);
    
    _pCylinderVector.push_back(_pLastCylinder);
    
    _pLastCylinder = c;
    
}

void Filament::PolymerizeBack(vector<double> coordinates) {
    
    Bead* b = BeadDB::Instance(BeadDBKey())->CreateBead(coordinates);
    Cylinder* c = CylinderDB::Instance(CylinderDBKey())->CreateCylinder(this, b);
    c->getMCylinder()->SetSecondBead(_pCylinderVector[0]->getMCylinder()->GetFirstBead());
    c->SetLast(false);
    _pCylinderVector.push_front(c);

}


void Filament::PolymerizeFront() {
    if (_pCylinderVector.size()<2) {cout<<"ERROR FILAMENT TO SHORT"<<endl;}
    
    else{
        
        auto tau = TwoPointDirection(
        _pCylinderVector[_pCylinderVector.size()-2]->getMCylinder()->GetFirstBead()->coordinate,
        _pCylinderVector[_pCylinderVector.size()-2]->getMCylinder()->GetSecondBead()->coordinate);
        
        Bead* b = BeadDB::Instance(BeadDBKey())->CreateBead(
                 NextPointProjection(_pCylinderVector[_pCylinderVector.size()-1]->getMCylinder()->GetSecondBead()->coordinate, L, tau) );
        Cylinder* c = CylinderDB::Instance(CylinderDBKey())->CreateCylinder(this, b, true);
        _pLastCylinder->getMCylinder()->SetSecondBead(b);
        _pLastCylinder->SetLast(false);
        _pCylinderVector.push_back(_pLastCylinder);
        _pLastCylinder = c;
                                       
        }


}

void Filament::PolymerizeBack() {
    if (_pCylinderVector.size()<2) {cout<<"ERROR FILAMENT TO SHORT"<<endl;}
    
    else{
        
        auto tau = TwoPointDirection(
         _pCylinderVector[0]->getMCylinder()->GetFirstBead()->coordinate, _pCylinderVector[0]->getMCylinder()->GetFirstBead()->coordinate);
        
        Bead* b = BeadDB::Instance(BeadDBKey())->CreateBead(
                    NextPointProjection(_pCylinderVector[0]->getMCylinder()->GetFirstBead()->coordinate, L, tau) );
        Cylinder* c = CylinderDB::Instance(CylinderDBKey())->CreateCylinder(this, b, true);
        c->getMCylinder()->SetSecondBead(_pCylinderVector[0]->getMCylinder()->GetFirstBead());
        c->SetLast(false);
        _pCylinderVector.push_front(c);
        
    }
    
    
}

vector<vector<double> > Filament::StraightFilamentProjection(vector<vector<double>> v, int numBeads){
    
    vector<vector<double>> coordinate;
    vector<double> tmpVec (3, 0);
    vector<double> tau (3, 0);
    double invD = 1/TwoPointDistance(v[1], v[0]);
    tau[0] = invD * ( v[1][0] - v[0][0] );
    tau[1] = invD * ( v[1][1] - v[0][1] );
    tau[2] = invD * ( v[1][2] - v[0][2] );
    
    for (int i = 0; i<numBeads; i++) {
        
        tmpVec[0] = v[0][0] + L * i * tau[0];
        tmpVec[1] = v[0][1] + L * i * tau[1];
        tmpVec[2] = v[0][2] + L * i * tau[2];
        
        coordinate.push_back(tmpVec);
    }
    return coordinate;
}


void Filament::DeleteBead(Bead*){
    
    cout<<"not implemented"<<endl;
}
