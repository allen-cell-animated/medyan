
//------------------------------------------------------------------
//  **M3SYM** - Simulation Package for the Mechanochemical
//              Dynamics of Active Networks, 3rd Generation
//
//  Copyright (2014) Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the Papoian lab page for more information:
//  http://papoian.chem.umd.edu/
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

#include "CompartmentGrid.h"

#include "MathFunctions.h"

using namespace mathfunc;

void BasicSnapshot::print(int step) {
    
    _outputFile.precision(10);
    
    // print first line (step number, time, number of filaments,
    // linkers, motors, branchers)
    _outputFile << step << " " << tau() << " " <<
        Filament::numFilaments() << " " <<
        Linker::numLinkers() << " " <<
        MotorGhost::numMotorGhosts() << " " <<
        BranchingPoint::numBranchingPoints() << endl;
    
    for(auto &filament : Filament::getFilaments()) {
        
        //print first line(Filament ID, length, left_delta, right_delta
        _outputFile << "F " << filament->getID() << " " <<
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
        _outputFile << "L " << linker->getID()<< " " <<
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
        _outputFile << "M " << motor->getID() << " " <<
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
        _outputFile << "B " << branch->getID() << " " <<
                               branch->getType() << endl;
        
        //print coordinates
        auto x = branch->coordinate;
        _outputFile<<x[0]<<" "<<x[1]<<" "<<x[2] << endl;
    }
    
    _outputFile <<endl;
}

void BirthTimes::print(int step) {
    
    _outputFile.precision(10);
    
    // print first line (step number, time, number of filaments,
    // linkers, motors, branchers)
    _outputFile << step << " " << tau() << " " <<
        Filament::numFilaments() << " " <<
        Linker::numLinkers() << " " <<
        MotorGhost::numMotorGhosts() << " " <<
        BranchingPoint::numBranchingPoints() << endl;
    
    for(auto &filament : Filament::getFilaments()) {
        
        //print first line(Filament ID, length, left_delta, right_delta
        _outputFile << "F " << filament->getID() << " " <<
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
        _outputFile << "L " << linker->getID()<< " " <<
                               linker->getType() << endl;
        
        //print birth times
        _outputFile << linker->getBirthTime() << " " <<
                       linker->getBirthTime() << endl;
    }
    
    for(auto &motor : MotorGhost::getMotorGhosts()) {
        
        //print first line
        _outputFile << "M " << motor->getID() << " " <<
                               motor->getType() << endl;
        
        //print birth times
        _outputFile << motor->getBirthTime() << " " <<
                       motor->getBirthTime() << endl;
    }
    
    for(auto &branch : BranchingPoint::getBranchingPoints()) {
        
        //print first line
        _outputFile << "B " << branch->getID() << " " <<
                               branch->getType() << endl;
        
        //print birth times
        _outputFile << branch->getBirthTime() << endl;
    }
    
    _outputFile <<endl;
}

void Forces::print(int step) {
    
    _outputFile.precision(10);
    
    // print first line (step number, time, number of filaments,
    // linkers, motors, branchers)
    _outputFile << step << " " << tau() << " " <<
        Filament::numFilaments() << " " <<
        Linker::numLinkers() << " " <<
        MotorGhost::numMotorGhosts() << " " <<
        BranchingPoint::numBranchingPoints() << endl;
    
    for(auto &filament : Filament::getFilaments()) {
        
        //print first line(Filament ID, length, left_delta, right_delta
        _outputFile << "F " << filament->getID() << " " <<
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
        _outputFile << "L " << linker->getID()<< " " <<
                               linker->getType() << endl;
        
        //print stretch force
        _outputFile << linker->getMLinker()->stretchForce << " " <<
                       linker->getMLinker()->stretchForce << endl;
    }
    
    for(auto &motor : MotorGhost::getMotorGhosts()) {
        
        //print first line
        _outputFile << "M " << motor->getID() << " " <<
                               motor->getType() << endl;
        
        //print stretch force
        _outputFile << motor->getMMotorGhost()->stretchForce << " " <<
                       motor->getMMotorGhost()->stretchForce << endl;
    }
    
    for(auto &branch : BranchingPoint::getBranchingPoints()) {
        
        //print first line
        _outputFile << "B " << branch->getID() << " " <<
                               branch->getType() << endl;
        
        //Nothing for branchers
        _outputFile << 0.0 << endl;
    }
    
    _outputFile <<endl;
}


void Stresses::print(int step) {

    _outputFile.precision(10);
    
    // print first line (step number, time, number of filaments,
    // linkers, motors, branchers)
    _outputFile << step << " " << tau() << " " <<
        Filament::numFilaments() << " " <<
        Linker::numLinkers() << " " <<
        MotorGhost::numMotorGhosts() << " " <<
        BranchingPoint::numBranchingPoints() << endl;
    
    for(auto &filament : Filament::getFilaments()) {
        
        //print first line(Filament ID, length, left_delta, right_delta
        _outputFile << "F " << filament->getID() << " " <<
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
        _outputFile << "L " << linker->getID()<< " " <<
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
        _outputFile << "M " << motor->getID() << " " <<
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
        _outputFile << "B " << branch->getID() << " " <<
                               branch->getType() << endl;
        
        //Nothing for branchers
        _outputFile << 0.0 << endl;
    }
    
    _outputFile <<endl;
}

void Chemistry::print(int step) {
    
    // print first line (step number, time)
    _outputFile << step << " " << tau() << endl;
    
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
    
    for(auto sf : _chemData.speciesFilament) {
        
        auto copyNum = Filament::countSpecies(sf);
        _outputFile << sf << ":FILAMENT " << copyNum << endl;
    }
    
    for(auto sp : _chemData.speciesPlusEnd) {
        
        auto copyNum = Filament::countSpecies(sp);
        _outputFile << sp << ":PLUSEND " << copyNum << endl;
    }
    
    for(auto sm : _chemData.speciesMinusEnd) {
        
        auto copyNum = Filament::countSpecies(sm);
        _outputFile << sm << ":MINUSEND " << copyNum << endl;
    }
    
    for(auto sl : _chemData.speciesLinker) {
        
        auto copyNum = Linker::countSpecies(sl);
        _outputFile << sl << ":LINKER " << copyNum << endl;
    }

    for(auto sm : _chemData.speciesMotor) {
        
        auto copyNum = MotorGhost::countSpecies(sm);
        _outputFile << sm << ":MOTOR " << copyNum << endl;
    }
    
    for(auto sb : _chemData.speciesBrancher) {
        
        auto copyNum = MotorGhost::countSpecies(sb);
        _outputFile << sb << ":BRANCHER " << copyNum << endl;
    }
    
    _outputFile <<endl;
}

