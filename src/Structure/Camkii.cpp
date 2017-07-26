#ifdef CAMKII
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



#include "SubSystem.h"
#include "CController.h"
#include "ChemManager.h"
#include "ChemRNode.h"

#include "Cylinder.h"
#include "Filament.h"
#include "Bead.h"

#include "GController.h"
#include "MathFunctions.h"
#include "Camkii.h"

using namespace mathfunc;

Database<Camkii*> Camkii::_camkiiDB; ///< Collection in SubSystem
void Camkii::updateCoordinate() {
    auto itr1 = _cylinders.begin();
    auto itr2 = _cylinders.begin() + 3;
    vector<double> avgMidpoint = {0.0, 0.0, 0.0};
    for(auto i=0; i<_cylinders.size()/2; i++){
        auto mp = midPointCoordinate(itr1[i]->coordinate, itr2[i]->coordinate, 0.5);
        avgMidpoint[0] += mp[0];
        avgMidpoint[1] += mp[1];
        avgMidpoint[2] += mp[2];
    }

    coordinate[0] = avgMidpoint[0]/3.0;
    coordinate[1] = avgMidpoint[1]/3.0;
    coordinate[2] = avgMidpoint[2]/3.0;
}

Camkii::Camkii(SubSystem* subsystem, int position, bool initialization)
    : _subSystem(subsystem), Trackable(true, true, true, false), _cylinders(), _position(position) {
    // TODO init cylinders
        
    vector<vector<double>> tmpBeadsCoord;
        
    Bead* b1 = _subSystem->addTrackable<Bead>(tmpBeadsCoord[0], this, 0);
    Bead* b2 = _subSystem->addTrackable<Bead>(tmpBeadsCoord[1], this, 1);
        
    Cylinder* c0 = _subSystem->addTrackable<Cylinder>(this, b1, b2, _filType, 0,
                                                          false, false, true);
        
    c0->setPlusEnd(true);
    c0->setMinusEnd(true);
    _cylinderVector.push_back(c0);
        
    for (int i = 2; i<numBeads; i++)
        extendPlusEnd(tmpBeadsCoord[i]);
        
        
    _ID = _camkiiDB.getID();
        
    if (initialization)
        cout<<"dddd";
    // TODO do we need this?
    parent->addChild(unique_ptr<Component>(this));
          
    //Set coordinate
    updateCoordinate();

    try {_compartment = GController::getCompartment(coordinate);}
    catch (exception& e) {
        cout << e.what() << endl;
        
        printSelf();
        
        exit(EXIT_FAILURE);
    }
        
    _compartment->addCamkii(this);
        
}

vector<vector<double>> Filament::zigZagFilamentProjection(vector<vector<double>>& v, int numBeads){
    
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
            tmpVec[0] = v[0][0] + SysParams::Geometry().cylinderSize[_filType] * i * tau[0];
            tmpVec[1] = v[0][1] + SysParams::Geometry().cylinderSize[_filType] * i * tau[1];
            tmpVec[2] = v[0][2] + SysParams::Geometry().cylinderSize[_filType] * i * tau[2];
        }
        else {
            tmpVec[0] = v[0][0] + SysParams::Geometry().cylinderSize[_filType] * i * perptau[0];
            tmpVec[1] = v[0][1] + SysParams::Geometry().cylinderSize[_filType] * i * perptau[1];
            tmpVec[2] = v[0][2] + SysParams::Geometry().cylinderSize[_filType] * i * perptau[2];
        }
        
        coordinate.push_back(tmpVec);
    }
    return coordinate;
}

Camkii::~Camkii() {
    //remove from compartment
    _compartment->removeCamkii(this);
    // TODO verify with James we need this
    for (auto &c: _cylinders){
        c->getCompartment()->removeCylinder(c);
    }
}

// TODO verify that's all we need with James.
void Camkii::updatePosition() {

    //check if were still in same compartment, set new position
    updateCoordinate();

    // TODO check if I need to update the binding manager
//
//    Compartment* c;
//
//    try {c = GController::getCompartment(coordinate);}
//    catch (exception& e) {
//        cout << e.what();
//
//        printSelf();
//
//        exit(EXIT_FAILURE);
//    }
//
//    if(c != _compartment) {
//
//#ifdef CHEMISTRY
//        auto oldCompartment = _compartment;
//        auto newCompartment = c;
//#endif
//
//        //remove from old compartment, add to new
//        _compartment->removeCamkii(this);
//        _compartment = c;
//        _compartment->addCamkii(this);
//
//#ifdef CHEMISTRY
//        auto oldCCylinder = _cCamkii.get();
//
//        //TODO camkii
//        //Remove old ccylinder from binding managers
//        for(auto &manager : oldCompartment->getFilamentBindingManagers())
//            manager->removePossibleBindings(oldCCylinder);
//
//        //clone and set new ccylinder
//        CCylinder* clone = _cCylinder->clone(c);
//        setCCylinder(clone);
//
//        auto newCCylinder = _cCylinder.get();
//
//        //Add new ccylinder to binding managers
//        for(auto &manager : newCompartment->getFilamentBindingManagers())
//            manager->addPossibleBindings(newCCylinder);
//    }
//#endif
//
//#ifdef MECHANICS
//    //update length
//    // TODO what to do here
//    _mCylinder->setLength(twoPointDistance(_b1->coordinate,
//                                           _b2->coordinate));
//#endif

}

/// @note -  The function uses the bead load force to calculate this changed rate.
/// If there is no force on the beads the reaction rates are set to the bare.

void Camkii::updateReactionRates() {
    // TODO double check if should be left empty
}

void Camkii::printSelf() {
    
    cout << endl;
    
    cout << "Camkii: ptr = " << this << endl;
    cout << "Camkii ID = " << _ID << endl;
    cout << "Parent ptr = " << getParent() << endl;
    cout << "Coordinates = " << coordinate[0] << ", " << coordinate[1] << ", " << coordinate[2] << endl;
    
    cout << "Position = " << _position << endl;
    cout << endl;

    cout << "Print cylinders..." << endl;
    for (const auto& c: _cylinders)
        c->printSelf();

    cout << endl;
}

// TODO
bool Camkii::isConsistent() {
    return true;
}

// TODO ????
// ChemManager* Camkii::_chemManager = 0; TODO what about that


#endif //CAMKII
