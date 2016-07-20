
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v3.0
//
//  Copyright (2015)  Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the MEDYAN web page for more information:
//  http://www.medyan.org
//------------------------------------------------------------------
#include <cmath>
#include <algorithm>

#include "Output.h"

#include "Filament.h"
#include "Cylinder.h"
#include "Bead.h"
#include "Linker.h"
#include "MotorGhost.h"
#include "BranchingPoint.h"
#include "Bubble.h"
#include "BoundaryElement.h" //added by jl135

#include "CompartmentGrid.h"

#include "MathFunctions.h"

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
        _outputFile << "FILAMENT " << filament->getID() << " " <<
        filament->getType() << " " <<
        filament->getCylinderVector().size() + 1 << " " <<
        filament->getDeltaMinusEnd() << " " << filament->getDeltaPlusEnd() << endl;

        //print coordinates
        for (auto cylinder : filament->getCylinderVector()){
            
            auto x = cylinder->getFirstBead()->coordinate;
            _outputFile<<x[0]<<" "<<x[1]<<" "<<x[2]<<" ";
            
        }
        //print last bead coord
        auto x = filament->getCylinderVector().back()->getSecondBead()->coordinate;
        _outputFile<<x[0]<<" "<<x[1]<<" "<<x[2];
        
        _outputFile << endl;
        
        //Reset deltas for this filament
        filament->resetDeltaPlusEnd();
        filament->resetDeltaMinusEnd();
    }
    
    
    for(auto &linker : Linker::getLinkers()) {
        
        //print first line
        _outputFile << "LINKER " << linker->getID()<< " " <<
                               linker->getType() << endl;
        
        //print coordinates
        auto x =
            midPointCoordinate(linker->getFirstCylinder()->getFirstBead()->coordinate,
                               linker->getFirstCylinder()->getSecondBead()->coordinate,
                               linker->getFirstPosition());
        _outputFile<<x[0]<<" "<<x[1]<<" "<<x[2] << " ";
        
        x = midPointCoordinate(linker->getSecondCylinder()->getFirstBead()->coordinate,
                               linker->getSecondCylinder()->getSecondBead()->coordinate,
                               linker->getSecondPosition());
        _outputFile<<x[0]<<" "<<x[1]<<" "<<x[2];
        
        _outputFile << endl;
    }

    for(auto &motor : MotorGhost::getMotorGhosts()) {
        
        //print first line
        _outputFile << "MOTOR " << motor->getID() << " " <<
                               motor->getType() << endl;
        
        //print coordinates
        auto x =
            midPointCoordinate(motor->getFirstCylinder()->getFirstBead()->coordinate,
                               motor->getFirstCylinder()->getSecondBead()->coordinate,
                               motor->getFirstPosition());
        _outputFile<<x[0]<<" "<<x[1]<<" "<<x[2] << " ";
        
        x = midPointCoordinate(motor->getSecondCylinder()->getFirstBead()->coordinate,
                               motor->getSecondCylinder()->getSecondBead()->coordinate,
                               motor->getSecondPosition());
        _outputFile<<x[0]<<" "<<x[1]<<" "<<x[2];
        
        _outputFile << endl;
    }
    
    for(auto &branch : BranchingPoint::getBranchingPoints()) {
        
        //print first line
        _outputFile << "BRANCHER " << branch->getID() << " " <<
                                      branch->getType() << endl;
        
        //print coordinates
        auto x = branch->coordinate;
        _outputFile<<x[0]<<" "<<x[1]<<" "<<x[2] << endl;
    }
    
    for(auto &bubble : Bubble::getBubbles()) {
        
        //print first line
        _outputFile << "BUBBLE " << bubble->getID() << " " <<
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
        _outputFile << "MOTOR " << motor->getID() << " " <<
                               motor->getType() << endl;
        
        //print birth times
        _outputFile << motor->getBirthTime() << " " <<
                       motor->getBirthTime() << endl;
    }
    
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
    // linkers, motors, branchers)
    _outputFile << snapshot << " " << tau() << " " <<
        Filament::numFilaments() << " " <<
        Linker::numLinkers() << " " <<
        MotorGhost::numMotorGhosts() << " " <<
        BranchingPoint::numBranchingPoints() << " " <<
        Bubble::numBubbles() << endl;
    
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
        _outputFile << "MOTOR " << motor->getID() << " " <<
                                   motor->getType() << endl;
        
        //print stretch force
        _outputFile << motor->getMMotorGhost()->stretchForce << " " <<
                       motor->getMMotorGhost()->stretchForce << endl;
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

    //added by jl135
    for (auto &be: BoundaryElement::getBoundaryElements()) {

    	//print first line
    	_outputFile << "TOTAL EXERTED FORCE BY ACTIN = " <<be->forceonboundaryAux<<endl;

    }
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
            
            _outputFile<< abs(k * deltaL) << " ";
            
        }
        //print last
        Cylinder* cylinder = filament->getCylinderVector().back();
        double k = cylinder->getMCylinder()->getStretchingConst();
        double deltaL = cylinder->getMCylinder()->getLength() -
                        cylinder->getMCylinder()->getEqLength();
        _outputFile<< abs(k * deltaL);
        
        _outputFile << endl;
    }
    
    for(auto &linker : Linker::getLinkers()) {
        
        //print first line
        _outputFile << "LINKER " << linker->getID()<< " " <<
                               linker->getType() << endl;
        
        //print
        double k = linker->getMLinker()->getStretchingConstant();
        double deltaL = linker->getMLinker()->getLength() -
                        linker->getMLinker()->getEqLength();
        
        
        _outputFile << abs(k * deltaL) << " " <<
                       abs(k * deltaL) << endl;
    }
    
    for(auto &motor : MotorGhost::getMotorGhosts()) {
        
        //print first line
        _outputFile << "MOTOR " << motor->getID() << " " <<
                               motor->getType() << endl;
        
        //print
        double k = motor->getMMotorGhost()->getStretchingConstant();
        double deltaL = motor->getMMotorGhost()->getLength() -
                        motor->getMMotorGhost()->getEqLength();
        
        _outputFile << abs(k * deltaL) << " " <<
                       abs(k * deltaL) << endl;
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
            
            auto copyNum = MotorGhost::countSpecies(sb);
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
/*//added by jl135, making a new class
void ForceOnBoundaryByActin::print(int snapshot) {

    _outputFile.precision(3);

    // print first line (step number, time)
    _outputFile << snapshot << " " << tau() << " " << endl;

    BoundaryElement::getBoundaryElements()->print(_outputFile);
    _outputFile << endl << endl;
}*/
