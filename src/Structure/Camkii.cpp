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



#include "CController.h"
#include "ChemManager.h"
#include "ChemRNode.h"

#include "Cylinder.h"
#include "Filament.h"
#include "Bead.h"

#include "GController.h"
#include "MathFunctions.h"
#include "Camkii.h"
#include "SysParams.h"

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

void Camkii::hexagonProjection(vector<double>& filamentInterfacePoint){

    vector<vector<double>> tmpBeadsCoord;
    double monoSize = SysParams::Geometry().monomerSize[_type];

    // pick random point which is within a certain distance away
    double directionX = Rand::randDouble(-1,1);
    double directionY = Rand::randDouble(-1,1);
    double directionZ = Rand::randDouble(-1,1);
    vector<double> randDirection1 = normalizedVector({directionX, directionY, directionZ});

    tmpBeadsCoord[0] = nextPointProjection(filamentInterfacePoint, monoSize/2, randDirection1);

    //switch rand direction and create second point
    randDirection1[0] = - randDirection1[0];
    randDirection1[1] = - randDirection1[1];
    randDirection1[2] = - randDirection1[2];

    tmpBeadsCoord[1] = nextPointProjection(tmpBeadsCoord[0], monoSize, randDirection1);

    directionX = Rand::randDouble(-1,1);
    directionY = Rand::randDouble(-1,1);
    directionZ = Rand::randDouble(-1,1);

    vector<double> randDirection2 = normalizedVector({directionX, directionY, directionZ});

    // sin(30) = sqrt(3)/2
    auto midPointParallelEdge = nextPointProjection(filamentInterfacePoint, sqrt(3)*monoSize, randDirection2);
    tmpBeadsCoord[4] = nextPointProjection(midPointParallelEdge, monoSize/2, randDirection1);

    //switch rand direction and create second point
    randDirection1[0] = - randDirection1[0];
    randDirection1[1] = - randDirection1[1];
    randDirection1[2] = - randDirection1[2];

    tmpBeadsCoord[4] = nextPointProjection(midPointParallelEdge, monoSize/2, randDirection1);
    tmpBeadsCoord[3] = nextPointProjection(tmpBeadsCoord[4], monoSize, randDirection1);

    tmpBeadsCoord[2] = nextPointProjection(midPointCoordinate(tmpBeadsCoord[3],tmpBeadsCoord[1], 0.5), monoSize/2,randDirection1);

    //switch rand direction and create second point
    randDirection1[0] = - randDirection1[0];
    randDirection1[1] = - randDirection1[1];
    randDirection1[2] = - randDirection1[2];

    tmpBeadsCoord[5] = nextPointProjection(midPointCoordinate(tmpBeadsCoord[4],tmpBeadsCoord[0], 0.5), monoSize/2,randDirection1);


    Bead* first = _subSystem->addTrackable<Bead>(tmpBeadsCoord[0], this, 0);
    Bead* curr = first;
    for (int i=0; i < tmpBeadsCoord.size() - 1; i++) {
        auto nextBeadCoord = tmpBeadsCoord[i+1];
        Bead* next = _subSystem->addTrackable<Bead>(nextBeadCoord, this, i+1);

        Cylinder* c = _subSystem->addTrackable<Cylinder>(this, curr, next, _type, 0,
                                                         false, false, true);
        c->setPlusEnd(false);
        c->setMinusEnd(false);
        _cylinders[i] = c;
        curr = next;
    }

    // Creating and adding the last cylinder
    Cylinder* c = _subSystem->addTrackable<Cylinder>(this, curr, first, _type, 0,
                                                     false, false, true);
    c->setPlusEnd(false);
    c->setMinusEnd(false);
    _cylinders[tmpBeadsCoord.size()-1] = c;

}

Camkii::Camkii(SubSystem* subsystem, short type, vector<double>& filamentInterfacePoint, bool initialization)
    : _subSystem(subsystem), _type(type),Trackable(true, true, true, false), _cylinders() {
    // TODO init cylinders
    hexagonProjection(filamentInterfacePoint);

    _ID = _camkiiDB.getID();
        
    if (initialization) // TODO what about initalization?
        cout<<"dddd";
    // TODO do we need this? do we need Composite?
    //parent->addChild(unique_ptr<Component>(this));
          
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

// TODO in the future.
// ChemManager* Camkii::_chemManager = 0; TODO what about that


#endif //CAMKII
