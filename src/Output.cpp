
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v3.2.1
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
#include "OutputStruct.h"

#include "Bead.h"
#include "BranchingPoint.h"
#include "Bubble.h"
#include "Cylinder.h"
#include "Filament.h"
#include "Linker.h"
#include "Membrane.hpp"
#include "MotorGhost.h"
#include "Vertex.h"

#include "Boundary.h"
#include "CompartmentGrid.h"
#include "core/controller/GController.h"
#include "Compartment.h"

#include "SysParams.h"
#include "MathFunctions.h"

using namespace mathfunc;

void BasicSnapshot::print(int snapshot) {
    
    _outputFile.precision(10);
    
    OutputStructSnapshot snapshots(snapshot);

    snapshots.outputFromSystem(_outputFile);
    
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
        _outputFile << "FILAMENT " << filament->getID() << " " <<
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
        _outputFile << "LINKER " << linker->getID()<< " " <<
                               linker->getType() << endl;
        
        //print birth times
        _outputFile << linker->getBirthTime() << " " <<
                       linker->getBirthTime() << endl;
    }
    
    for(auto &motor : MotorGhost::getMotorGhosts()) {
        
        //print first line
        //also contains a Bound(1) or unbound(0) qualifier
        _outputFile << "MOTOR " << motor->getID() << " " << motor->getType() << " " << 1 << endl;
        
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
        _outputFile << "BRANCHER " << branch->getID() << " " <<
                                      branch->getType() << endl;
        
        //print birth times
        _outputFile << branch->getBirthTime() << endl;
    }
    for(auto &bubble : Bubble::getBubbles()) {
        
        //print first line
        _outputFile << "BUBBLE " << bubble->getID() << " " <<
                                    bubble->getType() << endl;
        
        //print birth times
        _outputFile << bubble->getBead()->getBirthTime() << endl;
    }

    _outputFile <<endl;
}

void Forces::print(int snapshot) {
    
    _outputFile.precision(10);
    
    // print first line (snapshot number, time, number of filaments,
    // linkers, motors, branchers, membranes)
    _outputFile << snapshot << " " << tau() << " " <<
        Filament::numFilaments() << " " <<
        Linker::numLinkers() << " " <<
        MotorGhost::numMotorGhosts() << " " <<
        BranchingPoint::numBranchingPoints() << " " <<
        Bubble::numBubbles() << " " <<
        Membrane::numMembranes() << endl;
    
    for(auto &filament : Filament::getFilaments()) {
        
        //print first line (Filament ID, type, length, left_delta, right_delta
        _outputFile << "FILAMENT " << filament->getID() << " " <<
        filament->getType() << " " <<
        filament->getCylinderVector().size() + 1 << " " <<
        filament->getDeltaMinusEnd() << " " << filament->getDeltaPlusEnd() << endl;
        
        //print force
        for (auto cylinder : filament->getCylinderVector()){
            
            double forceMag= cylinder->getFirstBead()->FDotF();
            forceMag = sqrt(forceMag);
            _outputFile<<forceMag << " ";
            
        }
        //print last bead force
        double forceMag = filament->getCylinderVector().back()->
                          getSecondBead()->FDotF();
        forceMag = sqrt(forceMag);
        _outputFile<<forceMag;
        
        _outputFile << endl;
    }
    
    for(auto &linker : Linker::getLinkers()) {
        
        //print first line
        _outputFile << "LINKER " << linker->getID()<< " " <<
                               linker->getType() << endl;
        
        //print stretch force
        _outputFile << linker->getMLinker()->stretchForce << " " <<
                       linker->getMLinker()->stretchForce << endl;
    }
    
    for(auto &motor : MotorGhost::getMotorGhosts()) {
        
        //print first line
        //also contains a Bound(1) or unbound(0) qualifier
        _outputFile << "MOTOR " << motor->getID() << " " << motor->getType() << " " << 1 << endl;
        
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
        _outputFile << "BRANCHER " << branch->getID() << " " <<
                                      branch->getType() << endl;
        
        //Nothing for branchers
        _outputFile << 0.0 << endl;
    }
    for(auto &bubble : Bubble::getBubbles()) {
        
        //print first line
        _outputFile << "BUBBLE " << bubble->getID() << " " <<
                                    bubble->getType() << endl;
        
        //Nothing for bubbles
        _outputFile << 0.0 << endl;
    }
    
    for(auto &membrane : Membrane::getMembranes()) {
        
        //print first line (Membrane ID, type)
        _outputFile << "Membrane " << membrane->getId() << " " <<
            membrane->getType() << endl;
        
        //print force
        for(const auto& v : membrane->getMesh().getVertices()) {
            
            double forceMag = sqrt(v.attr.vertex->FDotF());
            _outputFile << forceMag << " ";
            
        }        
        _outputFile << endl;
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
        _outputFile << "FILAMENT " << filament->getID() << " " <<
        filament->getType() << " " <<
        filament->getCylinderVector().size() + 1 << " " <<
        filament->getDeltaMinusEnd() << " " << filament->getDeltaPlusEnd() << endl;
        
        //print
        for (auto cylinder : filament->getCylinderVector()){
            
            double k = cylinder->getMCylinder()->getStretchingConst();
            double deltaL = cylinder->getMCylinder()->getLength() -
                            cylinder->getMCylinder()->getEqLength();
            
            _outputFile<< k * deltaL << " ";
            
        }
        //print last
        Cylinder* cylinder = filament->getCylinderVector().back();
        double k = cylinder->getMCylinder()->getStretchingConst();
        double deltaL = cylinder->getMCylinder()->getLength() -
                        cylinder->getMCylinder()->getEqLength();
        _outputFile<< k * deltaL;
        
        _outputFile << endl;
    }
    
    for(auto &linker : Linker::getLinkers()) {
        
        //print first line
        _outputFile << "LINKER " << linker->getID()<< " " <<
                               linker->getType() << endl;
        
        //print
        _outputFile << linker->getMLinker()->stretchForce << " " <<
        linker->getMLinker()->stretchForce << endl;
    }
    
    for(auto &motor : MotorGhost::getMotorGhosts()) {
        
        //print first line
        //also contains a Bound(1) or unbound(0) qualifier
        _outputFile << "MOTOR " << motor->getID() << " " << motor->getType() << " " << 1 << endl;
        
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
        _outputFile << "BRANCHER " << branch->getID() << " " <<
                                      branch->getType() << endl;
        
        //Nothing for branchers
        _outputFile << 0.0 << endl;
    }
    for(auto &bubble : Bubble::getBubbles()) {
        
        //print first line
        _outputFile << "BUBBLE " << bubble->getID() << " " <<
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
        _outputFile << "FILAMENT " << filament->getID() << " " <<
        filament->getType() << " " <<
        filament->getCylinderVector().size() + 1 << " " <<
        filament->getDeltaMinusEnd() << " " << filament->getDeltaPlusEnd() << endl;
        
        //print
        for (auto cylinder : filament->getCylinderVector()){
            
            double k = SysParams::Mechanics().pinK;
            Bead* b = cylinder->getFirstBead();
            
            if(b->isPinned()) {
                auto norm = _subSystem->getBoundary()->normal(b->pinnedPosition);
                auto dirL = twoPointDirection(b->pinnedPosition, b->coordinate);
                
                double deltaL = twoPointDistance(b->coordinate, b->pinnedPosition);
                
                
                _outputFile<< k * deltaL * dotProduct(norm, dirL) << " ";
            }
            else
                _outputFile << 0.0 << " ";
            
        }
        //print last
        Cylinder* cylinder = filament->getCylinderVector().back();
        double k = SysParams::Mechanics().pinK;
        Bead* b = cylinder->getSecondBead();
        
        if(b->isPinned()) {
            auto norm = _subSystem->getBoundary()->normal(b->pinnedPosition);
            auto dirL = twoPointDirection(b->pinnedPosition, b->coordinate);
            
            double deltaL = twoPointDistance(b->coordinate, b->pinnedPosition);
            
            _outputFile<< k * deltaL * dotProduct(norm, dirL) << " ";
        }
        else
            _outputFile << 0.0 << " ";
        
        _outputFile << endl;
    }
    
    for(auto &linker : Linker::getLinkers()) {
        
        //print first line
        _outputFile << "LINKER " << linker->getID()<< " " <<
        linker->getType() << endl;
        
        _outputFile << 0.0 << " " << 0.0 << endl;
    }
    
    for(auto &motor : MotorGhost::getMotorGhosts()) {
        
        //print first line
        _outputFile << "MOTOR " << motor->getID() << " " <<
        motor->getType() << endl;
        
        _outputFile << 0.0 << " " << 0.0 << endl;
    }
    
    for(auto &branch : BranchingPoint::getBranchingPoints()) {
        
        //print first line
        _outputFile << "BRANCHER " << branch->getID() << " " <<
        branch->getType() << endl;
        
        //Nothing for branchers
        _outputFile << 0.0 << endl;
    }
    for(auto &bubble : Bubble::getBubbles()) {
        
        //print first line
        _outputFile << "BUBBLE " << bubble->getID() << " " <<
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
    Bubble::numBubbles() << endl;;
    
    for(auto &filament : Filament::getFilaments()) {
        
        //print first line (Filament ID, type, length, left_delta, right_delta)
        _outputFile << "FILAMENT " << filament->getID() << " " <<
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
        _outputFile << "LINKER " << linker->getID()<< " " <<
        linker->getType() << endl;
        
        _outputFile << linker->getType() << " " <<
        linker->getType() << endl;
    }
    
    for(auto &motor : MotorGhost::getMotorGhosts()) {
        
        //print first line
        //also contains a Bound(1) or unbound(0) qualifier
        _outputFile << "MOTOR " << motor->getID() << " " << motor->getType() << " " << 1 << endl;
        
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
        _outputFile << "BRANCHER " << branch->getID() << " " <<
        branch->getType() << endl;
        
        //Nothing for branchers
        _outputFile << branch->getType() << endl;
    }
    for(auto &bubble : Bubble::getBubbles()) {
        
        //print first line
        _outputFile << "BUBBLE " << bubble->getID() << " " <<
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
        _outputFile <<"FILAMENT " << filament->getID() << " " <<
        filament->getType() << " " <<
        filament->getCylinderVector().size() + 1 << " " <<
        filament->getDeltaMinusEnd() << " " << filament->getDeltaPlusEnd() << endl;

        //print plus end
        auto x = filament->getCylinderVector().back()->getSecondBead()->coordinate;
        _outputFile<<x[0]<<" "<<x[1]<<" "<<x[2]<<" \n";


        for (int i=0; i<filament->getCylinderVector().back()->getCCylinder()->getSize(); i++) {
            int out=filament->getCylinderVector().back()->getCCylinder()->getCMonomer(i)->activeSpeciesPlusEnd();
            if(out !=-1) {_outputFile << "PLUSEND: " << out << endl;}

        }

    }

    _outputFile << endl;

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
        _outputFile <<"FILAMENT " << filament->getID() << " " <<
        filament->getType() << " " <<
        filament->getCylinderVector().size() + 1 << " " <<
        filament->getDeltaMinusEnd() << " " << filament->getDeltaPlusEnd() << "\n"<<
        filament->getDeltaMinusEnd() << " " << filament->getDeltaPlusEnd() << " " <<
        filament->getPolyMinusEnd() << " " << filament->getPolyPlusEnd() << " " <<
        filament->getDepolyMinusEnd() << " " << filament->getDepolyPlusEnd() << " " <<
        filament->getNucleation() << " " << numMonomer << endl;

    }

    _outputFile << endl;

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

