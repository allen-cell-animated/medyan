
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v4.0
//
//  Copyright (2015-2018)  Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the MEDYAN web page for more information:
//  http://www.medyan.org
//------------------------------------------------------------------
#include <cmath>
#include <algorithm>

#include "Output.h"

#include "SubSystem.h"
#include "CompartmentGrid.h"

#include "Filament.h"
#include "Cylinder.h"
#include "Bead.h"
#include "Linker.h"
#include "MotorGhost.h"
#include "BranchingPoint.h"
#include "Bubble.h"

#include "Boundary.h"
#include "Compartment.h"
#include "GController.h"

#include "SysParams.h"
#include "MathFunctions.h"

#include "CController.h"
#include "ChemSimImpl.h"

using namespace mathfunc;

void BasicSnapshot::print(int snapshot) {

    _outputFile.precision(10);

    // print first line (snapshot number, time, number of filaments,
    // linkers, motors, branchers, bubbles)
    _outputFile << snapshot << " " << tau() << " " <<
        Filament::numFilaments() << " " <<
        Linker::numLinkers() << " " <<
        MotorGhost::numMotorGhosts() << " " <<
        BranchingPoint::numBranchingPoints() << " " <<
        Bubble::numBubbles() << endl;

    for(auto &filament : Filament::getFilaments()) {

        //print first line (Filament ID, type, length, left_delta, right_delta)
        _outputFile << "FILAMENT " << filament->getId() << " " <<
        filament->getType() << " " <<
        filament->getCylinderVector().size() + 1 << " " <<
        filament->getDeltaMinusEnd() << " " << filament->getDeltaPlusEnd() << endl;

        //print coordinates
        for (auto cylinder : filament->getCylinderVector()){

            auto x = cylinder->getFirstBead()->vcoordinate();
            _outputFile<<x[0]<<" "<<x[1]<<" "<<x[2]<<" ";

        }
        //print last bead coord
        auto x = filament->getCylinderVector().back()->getSecondBead()->vcoordinate();
        _outputFile<<x[0]<<" "<<x[1]<<" "<<x[2];

        _outputFile << endl;
    }


    for(auto &linker : Linker::getLinkers()) {

        //print first line
        _outputFile << "LINKER " << linker->getId()<< " " <<
                               linker->getType() << endl;

        //print coordinates
        auto x =
            midPointCoordinate(linker->getFirstCylinder()->getFirstBead()->vcoordinate(),
                               linker->getFirstCylinder()->getSecondBead()->vcoordinate(),
                               linker->getFirstPosition());
        _outputFile<<x[0]<<" "<<x[1]<<" "<<x[2] << " ";

        x = midPointCoordinate(linker->getSecondCylinder()->getFirstBead()->vcoordinate(),
                               linker->getSecondCylinder()->getSecondBead()->vcoordinate(),
                               linker->getSecondPosition());
        _outputFile<<x[0]<<" "<<x[1]<<" "<<x[2];

        _outputFile << endl;
    }

    for(auto &motor : MotorGhost::getMotorGhosts()) {

        //print first line
        //also contains a Bound(1) or unbound(0) qualifier
        _outputFile << "MOTOR " << motor->getId() << " " << motor->getType() << " " << 1 << endl;

        //print coordinates
        auto x =
            midPointCoordinate(motor->getFirstCylinder()->getFirstBead()->vcoordinate(),
                               motor->getFirstCylinder()->getSecondBead()->vcoordinate(),
                               motor->getFirstPosition());
        _outputFile<<x[0]<<" "<<x[1]<<" "<<x[2] << " ";

        x = midPointCoordinate(motor->getSecondCylinder()->getFirstBead()->vcoordinate(),
                               motor->getSecondCylinder()->getSecondBead()->vcoordinate(),
                               motor->getSecondPosition());
        _outputFile<<x[0]<<" "<<x[1]<<" "<<x[2];

        _outputFile << endl;
    }

    //DEPRECATED AS OF 9/8/16
//    //collect diffusing motors
//    for(auto md: _subSystem->getCompartmentGrid()->getDiffusingMotors()) {
//
//        int ID   = get<0>(md);
//        int type = get<1>(md);
//
//        auto firstPoint = get<2>(md);
//        auto secondPoint = get<3>(md);
//
//        _outputFile << "MOTOR " << ID << " " << type << " " << 0 << endl;
//
//        //print coordinates
//        _outputFile<<firstPoint[0]<<" "<<firstPoint[1]<<" "<<firstPoint[2] << " ";
//        _outputFile<<secondPoint[0]<<" "<<secondPoint[1]<<" "<<secondPoint[2];
//
//        _outputFile << endl;
//    }

    for(auto &branch : BranchingPoint::getBranchingPoints()) {

        //print first line
        _outputFile << "BRANCHER " << branch->getId() << " " <<
                                      branch->getType() << endl;

        //print coordinates
        auto x = branch->coordinate;
        _outputFile<<x[0]<<" "<<x[1]<<" "<<x[2] << endl;
    }

    for(auto &bubble : Bubble::getBubbles()) {

        //print first line
        _outputFile << "BUBBLE " << bubble->getId() << " " <<
                                    bubble->getType() << endl;

        //print coordinates
        auto x = bubble->coordinate;
        _outputFile<<x[0]<<" "<<x[1]<<" "<<x[2] << endl;
    }

    _outputFile <<endl;
}

void BirthTimes::print(int snapshot) {

    _outputFile.precision(10);

    // print first line (snapshot number, time, number of filaments,
    // linkers, motors, branchers, bubbles)
    _outputFile << snapshot << " " << tau() << " " <<
        Filament::numFilaments() << " " <<
        Linker::numLinkers() << " " <<
        MotorGhost::numMotorGhosts() << " " <<
        BranchingPoint::numBranchingPoints() << " " <<
        Bubble::numBubbles() << endl;

    for(auto &filament : Filament::getFilaments()) {

        //print first line (Filament ID, type, length, left_delta, right_delta)
        _outputFile << "FILAMENT " << filament->getId() << " " <<
        filament->getType() << " " <<
        filament->getCylinderVector().size() + 1 << " " <<
        filament->getDeltaMinusEnd() << " " << filament->getDeltaPlusEnd() << endl;

        //print birth times
        for (auto cylinder : filament->getCylinderVector()){

            auto b = cylinder->getFirstBead();
            _outputFile<< b->getBirthTime() << " ";

        }
        //last bead
        _outputFile<< filament->getCylinderVector().back()
                      ->getSecondBead()->getBirthTime();
        _outputFile << endl;
    }
    for(auto &linker : Linker::getLinkers()) {

        //print first line
        _outputFile << "LINKER " << linker->getId()<< " " <<
                               linker->getType() << endl;

        //print birth times
        _outputFile << linker->getBirthTime() << " " <<
                       linker->getBirthTime() << endl;
    }

    for(auto &motor : MotorGhost::getMotorGhosts()) {

        //print first line
        //also contains a Bound(1) or unbound(0) qualifier
        _outputFile << "MOTOR " << motor->getId() << " " << motor->getType() << " " << 1 << endl;
        //print birth times
        _outputFile << motor->getBirthTime() << " " <<
                       motor->getBirthTime() << endl;
    }

    //DEPRECATED AS OF 9/8/16
//
//    //collect diffusing motors
//    for(auto md: _subSystem->getCompartmentGrid()->getDiffusingMotors()) {
//
//        int ID   = get<0>(md);
//        int type = get<1>(md);
//
//        auto firstPoint = get<2>(md);
//        auto secondPoint = get<3>(md);
//
//        _outputFile << "MOTOR " << ID << " " << type << " " << 0 << endl;
//
//        //print coordinates
//        //print birth times
//        _outputFile << 0 << " " << 0 << endl;
//    }

    for(auto &branch : BranchingPoint::getBranchingPoints()) {

        //print first line
        _outputFile << "BRANCHER " << branch->getId() << " " <<
                                      branch->getType() << endl;

        //print birth times
        _outputFile << branch->getBirthTime() << endl;
    }
    for(auto &bubble : Bubble::getBubbles()) {

        //print first line
        _outputFile << "BUBBLE " << bubble->getId() << " " <<
                                    bubble->getType() << endl;

        //print birth times
        _outputFile << bubble->getBead()->getBirthTime() << endl;
    }

    _outputFile <<endl;
}

void Forces::print(int snapshot) {

    _outputFile.precision(10);

    // print first line (snapshot number, time, number of filaments,
    // linkers, motors, branchers)
    _outputFile << snapshot << " " << tau() << " " <<
        Filament::numFilaments() << " " <<
        Linker::numLinkers() << " " <<
        MotorGhost::numMotorGhosts() << " " <<
        BranchingPoint::numBranchingPoints() << " " <<
        Bubble::numBubbles() << endl;

    for(auto &filament : Filament::getFilaments()) {

        //print first line (Filament ID, type, length, left_delta, right_delta
        _outputFile << "FILAMENT " << filament->getId() << " " <<
        filament->getType() << " " <<
        filament->getCylinderVector().size() + 1 << " " <<
        filament->getDeltaMinusEnd() << " " << filament->getDeltaPlusEnd() << endl;

        //print force
        for (auto cylinder : filament->getCylinderVector()){
            
            floatingpoint forceMag= cylinder->getFirstBead()->FDotF();
            forceMag = sqrt(forceMag);
            _outputFile<<forceMag << " ";

        }
        //print last bead force
        floatingpoint forceMag = filament->getCylinderVector().back()->
                          getSecondBead()->FDotF();
        forceMag = sqrt(forceMag);
        _outputFile<<forceMag;

        _outputFile << endl;
    }

    for(auto &linker : Linker::getLinkers()) {

        //print first line
        _outputFile << "LINKER " << linker->getId()<< " " <<
                               linker->getType() << endl;

        //print stretch force
        _outputFile << linker->getMLinker()->stretchForce << " " <<
                       linker->getMLinker()->stretchForce << endl;
    }

    for(auto &motor : MotorGhost::getMotorGhosts()) {

        //print first line
        //also contains a Bound(1) or unbound(0) qualifier
        _outputFile << "MOTOR " << motor->getId() << " " << motor->getType() << " " << 1 << endl;
        //print stretch force
        _outputFile << motor->getMMotorGhost()->stretchForce << " " <<
                       motor->getMMotorGhost()->stretchForce << endl;
    }

    //DEPRECATED AS OF 9/8/16
//    //collect diffusing motors
//    for(auto md: _subSystem->getCompartmentGrid()->getDiffusingMotors()) {
//
//        int ID   = get<0>(md);
//        int type = get<1>(md);
//
//        auto firstPoint = get<2>(md);
//        auto secondPoint = get<3>(md);
//
//        _outputFile << "MOTOR " << ID << " " << type << " " << 0 << endl;
//
//        //print coordinates
//        //print birth times
//        _outputFile << 0 << " " << 0 << endl;
//    }

    for(auto &branch : BranchingPoint::getBranchingPoints()) {

        //print first line
        _outputFile << "BRANCHER " << branch->getId() << " " <<
                                      branch->getType() << endl;

        //Nothing for branchers
        _outputFile << 0.0 << endl;
    }
    for(auto &bubble : Bubble::getBubbles()) {

        //print first line
        _outputFile << "BUBBLE " << bubble->getId() << " " <<
                                    bubble->getType() << endl;

        //Nothing for bubbles
        _outputFile << 0.0 << endl;
    }

    _outputFile <<endl;
}

void Tensions::print(int snapshot) {

    _outputFile.precision(10);

    // print first line (snapshot number, time, number of filaments,
    // linkers, motors, branchers)
    _outputFile << snapshot << " " << tau() << " " <<
        Filament::numFilaments() << " " <<
        Linker::numLinkers() << " " <<
        MotorGhost::numMotorGhosts() << " " <<
        BranchingPoint::numBranchingPoints() << " " <<
        Bubble::numBubbles() << endl;;

    for(auto &filament : Filament::getFilaments()) {

        //print first line (Filament ID, type, length, left_delta, right_delta)
        _outputFile << "FILAMENT " << filament->getId() << " " <<
        filament->getType() << " " <<
        filament->getCylinderVector().size() + 1 << " " <<
        filament->getDeltaMinusEnd() << " " << filament->getDeltaPlusEnd() << endl;

        //print
        for (auto cylinder : filament->getCylinderVector()){
            
            floatingpoint k = cylinder->getMCylinder()->getStretchingConst();
            floatingpoint deltaL = cylinder->getMCylinder()->getLength() -
                            cylinder->getMCylinder()->getEqLength();

            _outputFile<< k * deltaL << " ";

        }
        //print last
        Cylinder* cylinder = filament->getCylinderVector().back();
        floatingpoint k = cylinder->getMCylinder()->getStretchingConst();
        floatingpoint deltaL = cylinder->getMCylinder()->getLength() -
                        cylinder->getMCylinder()->getEqLength();
        _outputFile<< k * deltaL;

        _outputFile << endl;
    }

    for(auto &linker : Linker::getLinkers()) {

        //print first line
        _outputFile << "LINKER " << linker->getId()<< " " <<
                               linker->getType() << endl;

        //print
        _outputFile << linker->getMLinker()->stretchForce << " " <<
        linker->getMLinker()->stretchForce << endl;
    }

    for(auto &motor : MotorGhost::getMotorGhosts()) {

        //print first line
        //also contains a Bound(1) or unbound(0) qualifier
        _outputFile << "MOTOR " << motor->getId() << " " << motor->getType() << " " << 1 << endl;
        //print
        _outputFile << motor->getMMotorGhost()->stretchForce << " " <<
        motor->getMMotorGhost()->stretchForce << endl;
    }

    //DEPRECATED AS OF 9/8/16
//    //collect diffusing motors
//    for(auto md: _subSystem->getCompartmentGrid()->getDiffusingMotors()) {
//
//        int ID   = get<0>(md);
//        int type = get<1>(md);
//
//        auto firstPoint = get<2>(md);
//        auto secondPoint = get<3>(md);
//
//        _outputFile << "MOTOR " << ID << " " << type << " " << 0 << endl;
//
//        //print coordinates
//        //print birth times
//        _outputFile << 0 << " " << 0 << endl;
//    }

    for(auto &branch : BranchingPoint::getBranchingPoints()) {

        //print first line
        _outputFile << "BRANCHER " << branch->getId() << " " <<
                                      branch->getType() << endl;

        //Nothing for branchers
        _outputFile << 0.0 << endl;
    }
    for(auto &bubble : Bubble::getBubbles()) {

        //print first line
        _outputFile << "BUBBLE " << bubble->getId() << " " <<
                                    bubble->getType() << endl;

        //Nothing for bubbles
        _outputFile << 0.0 << endl;
    }

    _outputFile <<endl;
}

void WallTensions::print(int snapshot) {

    _outputFile.precision(10);

    // print first line (snapshot number, time, number of filaments,
    // linkers, motors, branchers)
    _outputFile << snapshot << " " << tau() << " " <<
    Filament::numFilaments() << " " <<
    Linker::numLinkers() << " " <<
    MotorGhost::numMotorGhosts() << " " <<
    BranchingPoint::numBranchingPoints() << " " <<
    Bubble::numBubbles() << endl;;

    for(auto &filament : Filament::getFilaments()) {

        //print first line (Filament ID, type, length, left_delta, right_delta)
        _outputFile << "FILAMENT " << filament->getId() << " " <<
        filament->getType() << " " <<
        filament->getCylinderVector().size() + 1 << " " <<
        filament->getDeltaMinusEnd() << " " << filament->getDeltaPlusEnd() << endl;

        //print
        for (auto cylinder : filament->getCylinderVector()){
            
            floatingpoint k = SysParams::Mechanics().pinK;
            Bead* b = cylinder->getFirstBead();

            if(b->isPinned()) {
                auto norm = _subSystem->getBoundary()->normal(b->pinnedPosition);
                auto dirL = twoPointDirection(b->pinnedPosition, b->vcoordinate());
                
                floatingpoint deltaL = twoPointDistance(b->vcoordinate(), b->pinnedPosition);
                
                _outputFile<< k * deltaL * dotProduct(norm, dirL) << " ";
            }
            else
                _outputFile << 0.0 << " ";

        }
        //print last
        Cylinder* cylinder = filament->getCylinderVector().back();
        floatingpoint k = SysParams::Mechanics().pinK;
        Bead* b = cylinder->getSecondBead();

        if(b->isPinned()) {
            auto norm = _subSystem->getBoundary()->normal(b->pinnedPosition);
            auto dirL = twoPointDirection(b->pinnedPosition, b->vcoordinate());
            
            floatingpoint deltaL = twoPointDistance(b->vcoordinate(), b->pinnedPosition);
            
            _outputFile<< k * deltaL * dotProduct(norm, dirL) << " ";
        }
        else
            _outputFile << 0.0 << " ";

        _outputFile << endl;
    }

    for(auto &linker : Linker::getLinkers()) {

        //print first line
        _outputFile << "LINKER " << linker->getId()<< " " <<
        linker->getType() << endl;

        _outputFile << 0.0 << " " << 0.0 << endl;
    }

    for(auto &motor : MotorGhost::getMotorGhosts()) {

        //print first line
        _outputFile << "MOTOR " << motor->getId() << " " <<
        motor->getType() << endl;

        _outputFile << 0.0 << " " << 0.0 << endl;
    }

    for(auto &branch : BranchingPoint::getBranchingPoints()) {

        //print first line
        _outputFile << "BRANCHER " << branch->getId() << " " <<
        branch->getType() << endl;

        //Nothing for branchers
        _outputFile << 0.0 << endl;
    }
    for(auto &bubble : Bubble::getBubbles()) {

        //print first line
        _outputFile << "BUBBLE " << bubble->getId() << " " <<
        bubble->getType() << endl;

        //Nothing for bubbles
        _outputFile << 0.0 << endl;
    }

    _outputFile <<endl;
}

void Types::print(int snapshot) {

    _outputFile.precision(10);

    // print first line (snapshot number, time, number of filaments,
    // linkers, motors, branchers)
    _outputFile << snapshot << " " << tau() << " " <<
    Filament::numFilaments() << " " <<
    Linker::numLinkers() << " " <<
    MotorGhost::numMotorGhosts() << " " <<
    BranchingPoint::numBranchingPoints() << " " <<
    Bubble::numBubbles() << endl;

    for(auto &filament : Filament::getFilaments()) {

        //print first line (Filament ID, type, length, left_delta, right_delta)
        _outputFile << "FILAMENT " << filament->getId() << " " <<
        filament->getType() << " " <<
        filament->getCylinderVector().size() + 1 << " " <<
        filament->getDeltaMinusEnd() << " " << filament->getDeltaPlusEnd() << endl;

        //print
        for (auto cylinder : filament->getCylinderVector()){

            _outputFile<< cylinder->getType() << " ";

        }
        //print last
        Cylinder* cylinder = filament->getCylinderVector().back();
        _outputFile<< cylinder->getType();

        _outputFile << endl;
    }

    for(auto &linker : Linker::getLinkers()) {

        //print first line
        _outputFile << "LINKER " << linker->getId()<< " " <<
        linker->getType() << endl;

        _outputFile << linker->getType() << " " <<
        linker->getType() << endl;
    }

    for(auto &motor : MotorGhost::getMotorGhosts()) {

        //print first line
        //also contains a Bound(1) or unbound(0) qualifier
        _outputFile << "MOTOR " << motor->getId() << " " << motor->getType() << " " << 1 << endl;
        _outputFile << motor->getType() << " " <<
        motor->getType() << endl;
    }

    //DEPRECATED AS OF 9/8/16
//    //collect diffusing motors
//    for(auto md: _subSystem->getCompartmentGrid()->getDiffusingMotors()) {
//
//        int ID   = get<0>(md);
//        int type = get<1>(md);
//
//        auto firstPoint = get<2>(md);
//        auto secondPoint = get<3>(md);
//
//        _outputFile << "MOTOR " << ID << " " << type << " " << 0 << endl;
//
//        //print coordinates
//        //print birth times
//        _outputFile << type << " " << type << endl;
//    }

    for(auto &branch : BranchingPoint::getBranchingPoints()) {

        //print first line
        _outputFile << "BRANCHER " << branch->getId() << " " <<
        branch->getType() << endl;

        //Nothing for branchers
        _outputFile << branch->getType() << endl;
    }
    for(auto &bubble : Bubble::getBubbles()) {

        //print first line
        _outputFile << "BUBBLE " << bubble->getId() << " " <<
        bubble->getType() << endl;

        //Nothing for bubbles
        _outputFile << bubble->getType() << endl;
    }

    _outputFile <<endl;
}

void Chemistry::print(int snapshot) {

    // print first line (snapshot number, time)
    _outputFile << snapshot << " " << tau() << endl;

    // all diffusing and bulk species
    for(auto sd : _chemData.speciesDiffusing) {

        string name = get<0>(sd);
        auto copyNum = _grid->countDiffusingSpecies(name);

        _outputFile << name << ":DIFFUSING " << copyNum << endl;
    }

    for(auto sb : _chemData.speciesBulk) {

        string name = get<0>(sb);
        auto copyNum = _grid->countBulkSpecies(name);

        _outputFile << name << ":BULK " << copyNum << endl;
    }

    for(int filType = 0; filType < SysParams::Chemistry().numFilaments; filType++) {

        for(auto sf : _chemData.speciesFilament[filType]) {

            auto copyNum = Filament::countSpecies(filType, sf);
            _outputFile << sf << ":FILAMENT " << copyNum << endl;
        }

        for(auto sp : _chemData.speciesPlusEnd[filType]) {

            auto copyNum = Filament::countSpecies(filType, sp);
            _outputFile << sp << ":PLUSEND " << copyNum << endl;
        }

        for(auto sm : _chemData.speciesMinusEnd[filType]) {

            auto copyNum = Filament::countSpecies(filType, sm);
            _outputFile << sm << ":MINUSEND " << copyNum << endl;
        }

        for(auto sl : _chemData.speciesLinker[filType]) {

            auto copyNum = Linker::countSpecies(sl);
            _outputFile << sl << ":LINKER " << copyNum << endl;
        }

        for(auto sm : _chemData.speciesMotor[filType]) {

            auto copyNum = MotorGhost::countSpecies(sm);
            _outputFile << sm << ":MOTOR " << copyNum << endl;
        }

        for(auto sb : _chemData.speciesBrancher[filType]) {

            auto copyNum = BranchingPoint::countSpecies(sb);
            _outputFile << sb << ":BRANCHER " << copyNum << endl;
        }
    }

    _outputFile <<endl;
}

void MotorLifetimes::print(int snapshot) {

    _outputFile.precision(3);

    // print first line (snapshot number, time)
    _outputFile << snapshot << " " << tau() << " " << endl;

    MotorGhost::getLifetimes()->print(_outputFile);
    _outputFile << endl << endl;

    //clear list
    MotorGhost::getLifetimes()->clearValues();
}

void MotorWalkLengths::print(int snapshot) {

    _outputFile.precision(3);

    // print first line (snapshot number, time)
    _outputFile << snapshot << " " << tau() << " " << endl;

    MotorGhost::getWalkLengths()->print(_outputFile);
    _outputFile << endl << endl;

    //clear list
    MotorGhost::getWalkLengths()->clearValues();
}

void LinkerLifetimes::print(int snapshot) {

    _outputFile.precision(3);

    // print first line (step number, time)
    _outputFile << snapshot << " " << tau() << " " << endl;

    Linker::getLifetimes()->print(_outputFile);
    _outputFile << endl << endl;

    //clear list
    Linker::getLifetimes()->clearValues();
}

void FilamentTurnoverTimes::print(int snapshot) {

    _outputFile.precision(3);

    // print first line (step number, time)
    _outputFile << snapshot << " " << tau() << " " << endl;

    Filament::getTurnoverTimes()->print(_outputFile);
    _outputFile << endl << endl;
}

void Dissipation::print(int snapshot) {

    // print first line (snapshot number, time)
    _outputFile << snapshot << " " << tau() << endl;
    vector<floatingpoint> energies;
    energies = _cs->getEnergy();
    _outputFile << energies[0] << "     " << energies[1] << "     "<< energies[2]<<"     "<<energies[3]<<"     "<<energies[4];

    _outputFile <<endl;
}

void HRCD::print(int snapshot) {
    DissipationTracker* dt = _cs->getDT();
    vector<tuple<string, floatingpoint>> hrcdvec = dt->getHRCDVec();
    // print first line (snapshot number, time)

    _outputFile << snapshot << " " << tau() << endl;

    for(auto &i : hrcdvec){
        _outputFile<<get<0>(i)<<"     ";
    }
    _outputFile<<endl;
    for(auto &i : hrcdvec){
        _outputFile<<get<1>(i)<<"     ";
    }
    _outputFile<<endl<<endl;

}

void HRMD::print(int snapshot) {
    DissipationTracker* dt = _cs->getDT();
       // print first line (snapshot number, time)
    
    vector<tuple<string, floatingpoint>> cumHRMDMechEnergy = dt->getCumHRMDMechEnergy();
    vector<tuple<string, floatingpoint>> cumHRMDMechDiss = dt->getCumHRMDMechDiss();
    
    _outputFile << snapshot << " " << tau() << endl;
   
    // write row of names
    for(auto i = 0; i < cumHRMDMechEnergy.size(); i++){
        _outputFile << get<0>(cumHRMDMechEnergy[i]) << "     ";
    }
    _outputFile<<endl;
    
    // write row of mech energy
    for(auto i = 0; i < cumHRMDMechEnergy.size(); i++){
        _outputFile << get<1>(cumHRMDMechEnergy[i]) << "     ";
    }
    _outputFile<<endl;
    
    // write row of mech diss, assuming names are same
    for(auto i = 0; i < cumHRMDMechEnergy.size(); i++){
        for(auto j = 0; j < cumHRMDMechEnergy.size(); j++){
            if(get<0>(cumHRMDMechDiss[j]) == get<0>(cumHRMDMechEnergy[i])){
                _outputFile << get<1>(cumHRMDMechDiss[j]) << "     ";
            }
        }
        
    }
    
    _outputFile<<endl<<endl;
    
}

void PlusEnd::print(int snapshot) {

    _outputFile.precision(10);

    // print first line (snapshot number, time, number of filaments,
    // linkers, motors, branchers)
    _outputFile << snapshot << " " << tau() << " " <<
    Filament::numFilaments() << " " <<
    Linker::numLinkers() << " " <<
    MotorGhost::numMotorGhosts() << " " <<
    BranchingPoint::numBranchingPoints() << " " <<
    Bubble::numBubbles() <<endl;;

    for(auto &filament : Filament::getFilaments()) {

        //print first line (Filament ID, type, length, left_delta, right_delta)
        _outputFile <<"FILAMENT " << filament->getId() << " " <<
        filament->getType() << " " <<
        filament->getCylinderVector().size() + 1 << " " <<
        filament->getDeltaMinusEnd() << " " << filament->getDeltaPlusEnd() << endl;

        //print plus end
        auto x = filament->getCylinderVector().back()->getSecondBead()->vcoordinate();
        _outputFile<<x[0]<<" "<<x[1]<<" "<<x[2]<<" \n";


        for (int i=0; i<filament->getCylinderVector().back()->getCCylinder()->getSize(); i++) {
            int out=filament->getCylinderVector().back()->getCCylinder()->getCMonomer(i)->activeSpeciesPlusEnd();
            if(out !=-1) {_outputFile << "PLUSEND: " << out << endl;}

        }

    }

    _outputFile << endl;

}

void CMGraph::print(int snapshot) {

    // print first line (snapshot number, time)

    _outputFile << snapshot << " " << tau() << endl;

    vector<vector<int>> filIDVecL;
	vector<vector<int>> filIDVecM;
	vector<vector<int>> filIDVecB;

	map<uint64_t, array<int, 3>> filpaircounter;

	//Get filament pairs involved in each linker
    for(auto &linker : Linker::getLinkers()) {

        int fid1 = linker->getFirstCylinder()->getFilID();
        int fid2 = linker->getSecondCylinder()->getFilID();
	    vector<int> pair;
        if(fid1<fid2) {
	        pair.push_back(fid1);
	        pair.push_back(fid2);
        }
        else{
	        pair.push_back(fid2);
	        pair.push_back(fid1);
        }

        filIDVecL.push_back(pair);
    }

	//Get filament pairs involved in each motor
	for(auto &motor : MotorGhost::getMotorGhosts()) {

		int fid1 = motor->getFirstCylinder()->getFilID();
		int fid2 = motor->getSecondCylinder()->getFilID();
		vector<int> pair;
		if(fid1<fid2) {
			pair.push_back(fid1);
			pair.push_back(fid2);
		}
		else{
			pair.push_back(fid2);
			pair.push_back(fid1);
		}
		filIDVecM.push_back(pair);
	}

	//Get mother and daughter filament pairs involved in each brancher
	for(auto &brancher : BranchingPoint::getBranchingPoints()) {

		int fid1 = brancher->getFirstCylinder()->getFilID();
		int fid2 = brancher->getSecondCylinder()->getFilID();
		vector<int> pair;
		if(fid1<fid2) {
			pair.push_back(fid1);
			pair.push_back(fid2);
		}
		else{
			pair.push_back(fid2);
			pair.push_back(fid1);
		}
		filIDVecB.push_back(pair);
	}

	vector<vector<int>> uniqueFilIDVec;
    vector<array<int, 3>> uniqueFilIDVecCounts;

    for(auto i : filIDVecL){

        if(find(uniqueFilIDVec.begin(), uniqueFilIDVec.end(), i) != uniqueFilIDVec.end()) {

            int ind = find(uniqueFilIDVec.begin(), uniqueFilIDVec.end(), i) - uniqueFilIDVec.begin();

            uniqueFilIDVecCounts.at(ind)[0]++;

        } else {
            array<int, 3> pbVec={1,0,0};
            uniqueFilIDVecCounts.push_back(pbVec);
            uniqueFilIDVec.push_back(i);
        }
    }

    for(auto i : filIDVecM){

        if(find(uniqueFilIDVec.begin(), uniqueFilIDVec.end(), i) != uniqueFilIDVec.end()) {

            int ind = find(uniqueFilIDVec.begin(), uniqueFilIDVec.end(), i) - uniqueFilIDVec.begin();

            uniqueFilIDVecCounts.at(ind)[1]++;

        } else {
            array<int, 3> pbVec={0,1,0};
            uniqueFilIDVecCounts.push_back(pbVec);
            uniqueFilIDVec.push_back(i);
        }
    }

    for(auto i : filIDVecB){

        if(find(uniqueFilIDVec.begin(), uniqueFilIDVec.end(), i) != uniqueFilIDVec.end()) {

            int ind = find(uniqueFilIDVec.begin(), uniqueFilIDVec.end(), i) - uniqueFilIDVec.begin();

            uniqueFilIDVecCounts.at(ind)[2]++;

        } else {
            array<int, 3> pbVec={0,0,1};
            uniqueFilIDVecCounts.push_back(pbVec);
            uniqueFilIDVec.push_back(i);
        }
    }

    int count = 0;
    for(auto i: uniqueFilIDVecCounts){
        _outputFile<<uniqueFilIDVec.at(count)[0]<<" "<<uniqueFilIDVec.at(count)[1]<<" "<<
        i[0] <<" "<< i[1] << " " << i[2]<< " ";
        count++;
    }
    _outputFile<<endl<<endl;

}

void ReactionOut::print(int snapshot) {

    _outputFile.precision(10);

    // print first line (snapshot number, time, number of filaments,
    // linkers, motors, branchers)
    _outputFile << snapshot << " " << tau() << " " <<
    Filament::numFilaments() << " " <<
    Linker::numLinkers() << " " <<
    MotorGhost::numMotorGhosts() << " " <<
    BranchingPoint::numBranchingPoints() << " " <<
    Bubble::numBubbles() <<endl;;

    for(auto &filament : Filament::getFilaments()) {

        int numMonomer = 2; // 2 for plus/minus end
        for (auto c : filament->getCylinderVector()) {
            for (int i=0; i < c->getCCylinder()->getSize(); i++) {
                auto FilamentMonomer = c->getCCylinder()-> getCMonomer(i)->activeSpeciesFilament();
                if(FilamentMonomer != -1) {numMonomer ++;}

            }

        }

        //print first line (Filament ID, type, length, left_delta, right_delta)
        _outputFile <<"FILAMENT " << filament->getId() << " " <<
        filament->getType() << " " <<
        filament->getCylinderVector().size() + 1 << " " <<
        filament->getDeltaMinusEnd() << " " << filament->getDeltaPlusEnd() << "\n"<<
        filament->getDeltaMinusEnd() << " " << filament->getDeltaPlusEnd() << " " <<
        filament->getPolyMinusEnd() << " " << filament->getPolyPlusEnd() << " " <<
        filament->getDepolyMinusEnd() << " " << filament->getDepolyPlusEnd() << " " <<
        filament->getNucleation() << " " << numMonomer << endl;

        _outputFile << "SEVERING " << filament->getSevering() << endl;
        if (filament->getNewID().size() == 0) {
            _outputFile << "-1";
        }
        else {
            for (int i = 0; i < filament->getNewID().size(); ++i) {
                _outputFile << filament->getNewID()[i] << " ";
            }
        }

        _outputFile << endl;
    }

    _outputFile << endl;

}

void BRForces::print(int snapshot) {

    _outputFile.precision(10);

    // print first line (snapshot number, time, number of filaments,
    // linkers, motors, branchers)
    _outputFile << snapshot << " " << tau() << " " <<
    Filament::numFilaments() << " " <<
    Linker::numLinkers() << " " <<
    MotorGhost::numMotorGhosts() << " " <<
    BranchingPoint::numBranchingPoints() << " " <<
    Bubble::numBubbles() << endl;

    for(auto &filament : Filament::getFilaments()) {

        //print first line (Filament ID, type, length, left_delta, right_delta
        _outputFile << "FILAMENT " << filament->getId() << " " <<
        filament->getType() << " " <<
        filament->getCylinderVector().size() + 1 << " " <<
        filament->getDeltaMinusEnd() << " " << filament->getDeltaPlusEnd() << endl;

        //print force
        for (auto cylinder : filament->getCylinderVector()){

            floatingpoint forceMag= cylinder->getFirstBead()->brFDotbrF();
            forceMag = sqrt(forceMag);
            _outputFile<<forceMag << " ";

        }
        //print last bead force
        floatingpoint forceMag = filament->getCylinderVector().back()->
        getSecondBead()->brFDotbrF();
        forceMag = sqrt(forceMag);
        _outputFile<<forceMag;

        _outputFile << endl;
    }

    for(auto &linker : Linker::getLinkers()) {

        //print first line
        _outputFile << "LINKER " << linker->getId()<< " " <<
        linker->getType() << endl;

        //print stretch force
        _outputFile << linker->getMLinker()->stretchForce << " " <<
        linker->getMLinker()->stretchForce << endl;
    }

    for(auto &motor : MotorGhost::getMotorGhosts()) {

        //print first line
        //also contains a Bound(1) or unbound(0) qualifier
        _outputFile << "MOTOR " << motor->getId() << " " << motor->getType() << " " << 1 << endl;

        //print stretch force
        _outputFile << motor->getMMotorGhost()->stretchForce << " " <<
        motor->getMMotorGhost()->stretchForce << endl;
    }

}

void Concentrations::print(int snapshot) {

    _outputFile << snapshot << " " << tau() << endl;

    for(auto c : _subSystem->getCompartmentGrid()->getCompartments()) {

        if(c->isActivated()) {

            _outputFile << "COMPARTMENT: " << c->coordinates()[0] << " "
            << c->coordinates()[1] << " " << c->coordinates()[2] << endl;

            for(auto sd : _chemData.speciesDiffusing) {

                string name = get<0>(sd);
                auto s = c->findSpeciesByName(name);
                auto copyNum = s->getN();

                _outputFile << name << ":DIFFUSING " << copyNum << endl;
            }
        }
    }
    _outputFile << endl;
}

void MotorWalkingEvents::print(int snapshot) {
    DissipationTracker* dt = _cs->getDT();
    vector<tuple<floatingpoint, floatingpoint, floatingpoint, floatingpoint>> motorData = dt->getMotorData();
    for(auto i = 0; i < motorData.size(); i++){
        tuple<floatingpoint, floatingpoint, floatingpoint, floatingpoint> line = motorData[i];
        _outputFile<< get<0>(line) << "     " << get<1>(line) << "     "<< get<2>(line)<<"     "<<get<3>(line) <<endl;
    }
    dt->clearMotorData();


}

void LinkerUnbindingEvents::print(int snapshot) {
    DissipationTracker* dt = _cs->getDT();
    vector<tuple<floatingpoint, floatingpoint, floatingpoint, floatingpoint>> linkerUnbindingData = dt->getLinkerUnbindingData();
    for(auto i = 0; i < linkerUnbindingData.size(); i++){
        tuple<floatingpoint, floatingpoint, floatingpoint, floatingpoint> line = linkerUnbindingData[i];
        _outputFile<< get<0>(line) << "     " << get<1>(line) << "     "<< get<2>(line)<<"     "<<get<3>(line) <<endl;
    }
    dt->clearLinkerUnbindingData();


}

void LinkerBindingEvents::print(int snapshot) {
    DissipationTracker* dt = _cs->getDT();
    vector<tuple<floatingpoint, floatingpoint, floatingpoint, floatingpoint>> linkerBindingData = dt->getLinkerBindingData();
    for(auto i = 0; i < linkerBindingData.size(); i++){
        tuple<floatingpoint, floatingpoint, floatingpoint, floatingpoint> line = linkerBindingData[i];
        _outputFile<< get<0>(line) << "     " << get<1>(line) << "     "<< get<2>(line)<<"     "<<get<3>(line) <<endl;
    }
    dt->clearLinkerBindingData();


}

void Datadump::print(int snapshot) {
    _outputFile.close();
    _outputFile.open(_outputFileName, std::ofstream::trunc);
	_outputFile.precision(15);
    if(!_outputFile.is_open()) {
        cout << "There was an error opening file " << _outputFileName
             << " for output. Exiting." << endl;
        exit(EXIT_FAILURE);
    }
    //Rearrange bead and cylinder data to create a continuous array.
    Bead::rearrange();
	Cylinder::updateAllData();
    Cylinder::rearrange();
    _outputFile << snapshot << " " << tau() << endl;
    _outputFile << Filament::numFilaments() << " " <<
                             Linker::numLinkers() << " " <<
                             MotorGhost::numMotorGhosts() << " " <<
                             BranchingPoint::numBranchingPoints() << " " <<
                             Bubble::numBubbles() << endl;
    //Bead data
    _outputFile <<"BEAD DATA: BEADIDX(STABLE) FID COORDX COORDY COORDZ FORCEAUXX "
                  "FORCEAUXY FORCEAUXZ"<<endl;
    const auto& beadData = Bead::getDbDataConst();

    for(auto b:Bead::getBeads()){
        auto bidx = b->getStableIndex();
        Filament* f = static_cast<Filament*>(b->getParent());
        _outputFile <<bidx<<" "<<f->getId()<<" "<<beadData.coords.data()
        [3*bidx]<<" " <<beadData.coords.data()[3*bidx + 1]<<" "
        <<beadData.coords.data()[3*bidx + 2]<<" "<<beadData.forcesAux.data()[3*bidx]<<" "<<
        beadData.forcesAux.data()[3*bidx + 1]<<" "<<beadData.forcesAux.data()[3*bidx + 2]<<endl;

    }

    /*for(int bidx = 0; bidx<Bead::rawNumStableElements(); bidx++){

        _outputFile <<bidx<<" "<<beadData.coords.data()[3*bidx]<<" "<<beadData.coords.data()[3*bidx + 1]<<" "
        <<beadData.coords.data()[3*bidx + 2]<<" "<<beadData.forcesAux.data()[3*bidx]<<" "<<beadData
        .forcesAux.data()[3*bidx + 1]<<" "<<beadData.forcesAux.data()[3*bidx + 2]<<endl;
    }
	_outputFile <<endl;*/
    //Cylinder data
    _outputFile <<"CYLINDER DATA: CYLIDX(STABLE) B1_IDX B2_IDX MINUSENDTYPE "
                  "PLUSENDTYPE MINUSENDMONOMER PLUSENDMONOMER TOTALMONOMERS EQLEN"<<endl;
    const auto& cylinderInfoData = Cylinder::getDbData().value;

    for(int cidx = 0; cidx < Cylinder::rawNumStableElements(); cidx++){
        Cylinder* cyl = cylinderInfoData[cidx].chemCylinder->getCylinder();
        CCylinder* ccyl  = cylinderInfoData[cidx].chemCylinder;
        short filamentType = cyl->getType();
        int numMonomers = SysParams::Geometry().cylinderNumMon[filamentType];
        short minusendmonomer = 0;
        short plusendmonomer = numMonomers-1;
        short minusendtype = -1;
        short plusendtype = -1;
        short foundstatus = 0; //0 none found, 1 found one end, 2 found both ends
                for(int midx = 0; midx<numMonomers; midx++){
                if(foundstatus ==2)
                    break;
                short m = ccyl->getCMonomer(midx)->activeSpeciesMinusEnd();
                short p = ccyl->getCMonomer(midx)->activeSpeciesPlusEnd();
                if(m != -1) {
                    foundstatus++;
                    minusendtype = m;
                    minusendmonomer = midx;
                }

                if(p != -1) {
                    plusendtype = p;
                    foundstatus++;
                    plusendmonomer = midx;
                }
            }

        /*Cidx minus-end plus-end num-monomers*/
        _outputFile <<cidx<<" "<<cyl->getFirstBead()->getStableIndex()<<" "
        <<cyl->getSecondBead()->getStableIndex()<<" "<<minusendtype<<" "<<plusendtype<<" "
        <<minusendmonomer<<" "<<plusendmonomer<<" "<<(plusendmonomer-minusendmonomer)+1<<" "
        <<cyl->getMCylinder()->getEqLength()<<endl;
    }
	_outputFile <<endl;
    //Filament Data
	_outputFile <<"FILAMENT DATA: FILID CYLIDvec"<<endl;
	for(auto fil : Filament::getFilaments()){
		_outputFile <<fil->getId()<<" ";
		for(auto cyl :fil->getCylinderVector()){
			_outputFile << cyl->getStableIndex()<<" ";
		}
		_outputFile << endl;
	}
	_outputFile <<endl;
	//Linker Data
	_outputFile <<"LINKER DATA: LINKERID CYL1_IDX CYL2_IDX POS1 POS2 EQLEN"<<endl;
	for(auto l :Linker::getLinkers()){
		Cylinder* cyl1 = l->getFirstCylinder();
		Cylinder* cyl2 = l->getSecondCylinder();
		float pos1 = l->getFirstPosition();
		float pos2 = l->getSecondPosition();
		_outputFile <<l->getId()<<" "<<cyl1->getStableIndex()<<" "<<cyl2->getStableIndex
		()<<" "<<pos1<<" "<<pos2<<" "<<l->getMLinker()->getEqLength()<<endl;
	}
	_outputFile <<endl;
	//MOTOR Data
	_outputFile <<"MOTOR DATA: MOTORID CYL1_IDX CYL2_IDX POS1 POS2 EQLEN"<<endl;
	for(auto l :MotorGhost::getMotorGhosts()){
		Cylinder* cyl1 = l->getFirstCylinder();
		Cylinder* cyl2 = l->getSecondCylinder();
		float pos1 = l->getFirstPosition();
		float pos2 = l->getSecondPosition();
		_outputFile <<l->getId()<<" "<<cyl1->getStableIndex()<<" "<<cyl2->getStableIndex
				()<<" "<<pos1<<" "<<pos2<<" "<<l->getMMotorGhost()->getEqLength()<<endl;
	}
	_outputFile <<endl;
	//Brancher Data
	_outputFile <<"BRANCHING DATA: BRANCHID CYL1_IDX CYL2_IDX POS1 EQLEN"<<endl;
	for(auto l :BranchingPoint::getBranchingPoints()){
		Cylinder* cyl1 = l->getFirstCylinder();
		Cylinder* cyl2 = l->getSecondCylinder();
		float pos1 = l->getPosition();
		_outputFile <<l->getId()<<" "<<cyl1->getStableIndex()<<" "<<cyl2->getStableIndex
				()<<" "<<pos1<<" "<<l->getMBranchingPoint()->getEqLength()<<endl;
	}
	_outputFile <<endl;
	//Compartment Data
	_outputFile <<"COMPARTMENT DATA: CMPID DIFFUSINGSPECIES"
			   "COPYNUM"<<endl;
	for(auto cmp:_subSystem->getCompartmentGrid()->getCompartments()){
		_outputFile <<cmp->getId()<<" ";
		for(auto sd : _chemData.speciesDiffusing) {
			string name = get<0>(sd);
			auto s = cmp->findSpeciesByName(name);
			auto copyNum = s->getN();
			_outputFile <<name<<" "<<copyNum<<" ";
		}

		_outputFile <<endl;
	}
	_outputFile <<endl;
	//BulkSpecies
	_outputFile <<"BULKSPECIES: BULKSPECIES COPYNUM"<<endl;
	auto cmp = _subSystem->getCompartmentGrid()->getCompartments()[0];
	for(auto sb : _chemData.speciesBulk) {
		string name = get<0>(sb);
		auto s = cmp->findSpeciesByName(name);
		auto copyNum = s->getN();
		_outputFile <<name<<" "<<copyNum<<" ";
	}
	_outputFile <<endl;
}

void HessianMatrix::print(int snapshot){
    _outputFile.precision(10);
    vector<vector<vector<floatingpoint > > > hVec = _ffm->hessianVector;
    vector<floatingpoint> tauVector = _ffm-> tauVector;
    // Outputs a sparse representation of the Hessian matrix, where only elements with appreciable size (>0.00001) are
    // output along with their indices.  Currently this outputs for each minimization, however to reduce the file size this could be changed.
    for(auto k = 0; k < hVec.size(); k++){
        vector<vector<floatingpoint > > hMat = hVec[k];
        int total_DOF = hMat.size();
        vector<tuple<int, int, floatingpoint>> elements;
        for(auto i = 0; i < total_DOF; i++){
            for(auto j = 0; j < total_DOF; j++){
                if(std::abs(hMat[i][j]) > 0.00001){
                    elements.push_back(std::make_tuple(i,j,hMat[i][j]));
                }
            }
        }
        _outputFile << tauVector[k] << "     "<< total_DOF<< "     " << elements.size()<<endl;
        for(auto i = 0; i < elements.size(); i++){
            tuple<int, int, floatingpoint> element = elements[i];
            _outputFile<< get<0>(element) << "     "<< get<1>(element)<<"     "<< get<2>(element)<<endl;
        }
    _outputFile<<endl;
    };
    // This clears the vectors storing the matrices to reduce the amount of memory needed.  
    _ffm->clearHessian();
}

void TwoFilament::print(int snapshot){
    _outputFile.precision(10);

    vector<vector<floatingpoint>> coordfil1;
    vector<vector<floatingpoint>> coordfil2;
    int count = 0;

    for(auto &filament : Filament::getFilaments()) {
        for (auto cylinder : filament->getCylinderVector()){
            auto x = cylinder->getFirstBead()->vcoordinate();
            if(count == 0)
            	coordfil1.push_back(x);
            else
            	coordfil2.push_back(x);
        }
        count++;
    }
    if(coordfil1.size() != coordfil1.size()) {
        cout << "Sizes do not match! Exiting. Reinitialize with same number of beads in "
                "each filament" << endl;
        exit(EXIT_FAILURE);
    }
    auto Nbeads = coordfil1.size();
    double avgdissq = 0;
    for(auto i = 0; i < Nbeads; i ++){
    	auto c1 = coordfil1[Nbeads-1-i];
    	auto c2 = coordfil2[i];
    	avgdissq += twoPointDistancesquared(c1, c2);
    }
    _outputFile << tau()<< " "<<sqrt(avgdissq)<<endl;
}

