//
//  MFilament.cpp
//  CytoMech
//
//  Created by Konstantin Popov on 4/15/14.
//  Copyright (c) 2014 Konstantin Popov. All rights reserved.
//

#include "MFilament.h"

#include "MFilament.h"
#include "MBead.h"
#include "MSystem.h"

using namespace std;
using namespace mathfunc;



Filament::Filament(System* ps, Network* pn ,vector<double> position, vector<double> direction){
   
    
    //This constructor creates a short filament, containing only two beads. Coordinates of the first bead is an input, second is set up by using an input direction and a coarsegrain length L. Using all this, two constructors for beads are called.
    
    _pSystem = ps;
    _pNetwork = pn;
	
    
    Bead* b1 = BeadDB::Instance()->AddNewBead(position);
    Cylinder* c1 = CylinderDB::Instance()->AddNewCylinder(this, b1);
    _pLastCylinder = c1;
   

    
    PolymerizeFront( NextPointProjection(position, L, direction) );
}


Filament::Filament(System* ps, Network* pn, vector<vector<double> > position, int numBeads){
    
    /// This constructor is called to create a longer filament. It creates a filament with a number of beads numBeads. Filaments starts and ends in the point determined by position vector and has a direction direction. Number of beads is equal to the number of cylinders. The last cylinder doesnt have an end(second) bead and will not be pushed to cylinder vector, but will be stored in the _pLastCylinder;
    
    _pSystem = ps;
    _pNetwork = pn;
    
    vector<vector<double> > tmpBeadsCoord = StraightFilamentProjection(position, numBeads); //this function calculate coordinates for all beads on the line separated by a segment length.
   
    Bead* b0 = BeadDB::Instance(BeadDBKey())->AddNewBead(tmpBeadsCoord[0]);
    Cylinder* c0 = CylinderDB::Instance(CylinderDBKey())->AddNewCylinder(this, b0);
    _pLastCylinder = c0;
   
    
    for (int i = 1; i<numBeads; i++) {
        
        PolymerizeFront(tmpBeadsCoord[i]);  //Create n beads and n cylinders: x---x----x...x----x----o.
        }
    
}

void Filament::PolymerizeFront(vector<double> coordinates) {
    
    Bead* b = BeadDB::Instance(BeadDBKey())->AddNewBead(coordinates);
    Cylinder* c = CylinderDB::Instance(CylinderDBKey())->AddNewCylinder(this, b);
    _pLastCylinder->getMCylinder()->SetSecondBead(b);
    _pLastCylinder->SetLast(false);
    
    _pCylinderVector.push_back(_pLastCylinder);
    
    _pLastCylinder = c;
    
}

void Filament::PolymerizeBack(vector<double> coordinates) {
    
    Bead* b = BeadDB::Instance(BeadDBKey())->AddNewBead(coordinates);
    Cylinder* c = CylinderDB::Instance(CylinderDBKey())->AddNewCylinder(this, b);
    c->SetSecondBead(_pCylinderVector[0]->getMCylinder()->GetFirstBead());
    c->SetLast(false);
    _pCylinderVector.push_front(c);

}


void Filament::PolymerizeFront() {
    if (_pCylinderVector.size()<2) {cout<<"ERROR FILAMENT TO SHORT"<<endl;}
    
    else{
        
        auto tau = TwoPointDirection(_pCylinderVector[_pCylinderVector.size()-2]->getMCylinder()->GetFirstBead()->coordinate, _pCylinderVector[_pCylinderVector.size()-2]->getMCylinder()->GetSecondBead()->coordinate);
        
        Bead* b = BeadDB::Instance(BeadDBKey())->AddNewBead( NextPointProjection(_pCylinderVector[_pCylinderVector.size()-1]->getMCylinder()->GetSecondBead()->coordinate, L, tau) );
        Cylinder* c = CylinderDB::Instance(CylinderDBKey())->AddNewCylinder(this, b, false);
        _pLastCylinder->getMCylinder()->SetSecondBead(b);
        _pLastCylinder->SetLast(false);
        _pCylinderVector.push_back(_pLastCylinder);
        _pLastCylinder = c;
                                       
        }


}

void Filament::PolymerizeBack() {
    if (_pCylinderVector.size()<2) {cout<<"ERROR FILAMENT TO SHORT"<<endl;}
    
    else{
        
        auto tau = TwoPointDirection(_pCylinderVector[0]->getMCylinder()->GetFirstBead()->coordinate, _pCylinderVector[0]->getMCylinder()->GetFirstBead()->coordinate);
        
        Bead* b = BeadDB::Instance( NextPointProjection(_pCylinderVector[0]->getMCylinder()->GetFirstBead()->coordinate, L, tau) );
        Cylinder* c = CylinderDB::Instance(CylinderDBKey())->AddNewCylinder(this, b, false);
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

//Public mechanical interfaces which call private "potentials" and force expressions to calculate bonded interactions:

double Filament::CopmuteEnergy(){
    
    ///Iterates over beads in the filament and calculate energy from stretchin, bending and twisting interaction.
    
    double U = 0;
    double Us = 0;
    double Ub = 0;
    double Ut = 0;
    
    for ( auto it = _pCylinderVector.begin()+1; it != _pCylinderVector.end(); it++){
        
        Us += EnergyHarmonicStretching( (*(it-1)) );
        Ut += EnergyHarmonicBending( (*(it-1)), (*it) );
        
    }
    
    U = Us + Ub + Ut;
    
    return U;
    
}

double Filament::CopmuteEnergy(double d){
    
    ///Iterates over beads in the filament and calculate energy from stretchin, bending and twisting interaction(with coord - d*force argument).
    
    double U = 0;
    double Us = 0;
    double Ub = 0;
    double Ut = 0;
    
    for ( auto it = _pCylinderVector.begin()+1; it != _pCylinderVector.end(); it++){
        
        //        cout<<"force size  " <<(*(it-1))->force.size()<<endl;
        //        cout<<"coord size  " <<(*(it-1))->coordinate.size()<<endl;
        
        Us += EnergyHarmonicStretching((*(it-1)), d);
        Ut += EnergyHarmonicBending( (*(it-1)), (*it), d );
        
    }
    
    U = Us + Ub + Ut;
    
    //    cout<<"Uaux =  "<<U <<endl;
    return U;
}

/// Interface that iterates over beads in a filament and compute forces on aech pair(triplet) of beads i and i-1 (i-1, i, i+1; );
void Filament:: CopmuteForce(){
    
    cout<< "fz= "<<_pCylinderVector[0]->GetFirstBead()->force[2]<<endl;
    
    for ( auto it = _pCylinderVector.begin()+1; it != _pCylinderVector.end(); it++){
        
        ForceHarmonicStretching((*(it-1)));
        
        ForceHarmonicBending( (*(it-1)), (*it) );
        
    }
    
    cout<< "fz= "<<_pCylinderVector[0]->GetFirstBead()->force[2]<<endl;
    
}


void Filament:: CopmuteForceAux(){
    
   for ( auto it = _pCylinderVector.begin()+1; it != _pCylinderVector.end(); it++){
        
       ForceHarmonicStretchingAux((*(it-1)));
       
       ForceHarmonicBendingAux( (*(it-1)), (*it) );
       
        
    }
    
}

void Filament::DeleteBead(Bead*){
    
    cout<<"not implemented"<<endl;
}
