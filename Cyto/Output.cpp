//
//  Output.cpp
//  Cyto
//
//  Created by James Komianos on 9/16/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#include "Output.h"
#include "MathFunctions.h"

using namespace mathfunc;

///Print basic information about filaments
void Output::printBasicSnapshot(int step) {
    
    _outputFile.precision(10);
    
    //print first line (step number, time, number of filaments, linkers, motors)
    _outputFile << step << " " << tau() << " " << FilamentDB::Instance(FilamentDBKey())->size() <<
    " " << LinkerDB::Instance(LinkerDBKey())->size() << " " << MotorGhostDB::Instance(MotorGhostDBKey())->size() << std::endl;
    
    for(auto &filament : *FilamentDB::Instance(FilamentDBKey())) {

        ///print first line(Filament ID, length, left_delta, right_delta
        _outputFile << "F " << filament->getID() << " " << filament->getCylinderVector().size() + 1
            << " " << filament->getDeltaMinusEnd() << " " << filament->getDeltaPlusEnd() << std::endl;

        ///print coordinates
        for (auto cylinder : filament->getCylinderVector()){
            
            auto x = cylinder->GetFirstBead()->coordinate;
            _outputFile<<x[0]<<" "<<x[1]<<" "<<x[2]<<" ";
            
        }
        ///print last bead coord
        auto x = filament->getCylinderVector().back()->GetSecondBead()->coordinate;
        _outputFile<<x[0]<<" "<<x[1]<<" "<<x[2];
        
        _outputFile << std::endl;
        
        ///Reset deltas for this filament
        filament->resetDeltaPlusEnd();
        filament->resetDeltaMinusEnd();
    }
    
    
    for(auto &linker : *LinkerDB::Instance(LinkerDBKey())) {
        
        ///print first line(Filament ID, length, left_delta, right_delta
        _outputFile << "L " << linker->getLinkerID()<< " " << linker->getLinkerType() << std::endl;
        
        ///print coordinates
        auto x = MidPointCoordinate(linker->getFirstCylinder()->GetFirstBead()->coordinate,
                                    linker->getFirstCylinder()->GetSecondBead()->coordinate,
                                    linker->getFirstPosition());
        _outputFile<<x[0]<<" "<<x[1]<<" "<<x[2] << " ";
        
        x = MidPointCoordinate(linker->getSecondCylinder()->GetFirstBead()->coordinate,
                                    linker->getSecondCylinder()->GetSecondBead()->coordinate,
                                    linker->getSecondPosition());
        _outputFile<<x[0]<<" "<<x[1]<<" "<<x[2];
        
        _outputFile << std::endl;
    }

    for(auto &motor : *MotorGhostDB::Instance(MotorGhostDBKey())) {
        
        ///print first line(Filament ID, length, left_delta, right_delta
        _outputFile << "M " << motor->getMotorID() << " " << motor->getMotorType() << std::endl;
        
        ///print coordinates
        auto x = MidPointCoordinate(motor->getFirstCylinder()->GetFirstBead()->coordinate,
                                    motor->getFirstCylinder()->GetSecondBead()->coordinate,
                                    motor->getFirstPosition());
        _outputFile<<x[0]<<" "<<x[1]<<" "<<x[2] << " ";
        
        x = MidPointCoordinate(motor->getSecondCylinder()->GetFirstBead()->coordinate,
                               motor->getSecondCylinder()->GetSecondBead()->coordinate,
                               motor->getSecondPosition());
        _outputFile<<x[0]<<" "<<x[1]<<" "<<x[2];
        
        _outputFile << std::endl;
    }

    
    _outputFile <<std::endl;
}


void Output::printSnapshot(int step) {
    
    _outputFile.precision(10);
    
    //print first line (step number, time, number of filaments
    _outputFile << step << " " << tau() << " " << FilamentDB::Instance(FilamentDBKey())->size() << std::endl;
    
    for(auto &filament : *FilamentDB::Instance(FilamentDBKey())) {
        
        ///print first line(Filament ID, length, index of first bead, index of last bead
        _outputFile << filament->getID() << " " << filament->getCylinderVector().size() + 1
        << " " << filament->getCylinderVector().front()->GetFirstBead()->getID()
        << " " << filament->getCylinderVector().back()->GetSecondBead()->getID() << std::endl;
        
        ///print coordinates
        for (auto cylinder : filament->getCylinderVector()){
            
            auto x = cylinder->GetFirstBead()->coordinate;
            _outputFile<<x[0]<<" "<<x[1]<<" "<<x[2]<<" ";
            
        }
        auto x = filament->getCylinderVector().back()->GetSecondBead()->coordinate;
        _outputFile<<x[0]<<" "<<x[1]<<" "<<x[2];
        
        _outputFile << std::endl;
    }
    _outputFile <<std::endl;
    
}


void Output::printBirthTimes(int step) {
    
    _outputFile.precision(10);
    
    //print first line (step number, time, number of filaments
    _outputFile << step << " " << tau() << " " << FilamentDB::Instance(FilamentDBKey())->size() << std::endl;
    
    for(auto &filament : *FilamentDB::Instance(FilamentDBKey())) {
        
        ///print first line(Filament ID, length, index of first bead, index of last bead
        _outputFile << filament->getID() << " " << filament->getCylinderVector().size() + 1
        << " " << filament->getCylinderVector().front()->GetFirstBead()->getID()
        << " " << filament->getCylinderVector().back()->GetSecondBead()->getID() << std::endl;
        
        ///print coordinates
        for (auto cylinder : filament->getCylinderVector()){
            
            auto b = cylinder->GetFirstBead();
            _outputFile<< b->getBirthTime() << " ";
            
        }
        ///last bead
        _outputFile<< filament->getCylinderVector().back()->GetSecondBead()->getBirthTime();
        _outputFile << std::endl;
    }
    _outputFile <<std::endl;
}

void Output::printForces(int step) {
    
    _outputFile.precision(10);
    
    //print first line (step number, time, number of filaments
    _outputFile << step << " " << tau() << " " << FilamentDB::Instance(FilamentDBKey())->size() << std::endl;
    
    for(auto &filament : *FilamentDB::Instance(FilamentDBKey())) {
        
        ///print first line(Filament ID, length, index of first bead, index of last bead
        _outputFile << filament->getID() << " " << filament->getCylinderVector().size() + 1
        << " " << filament->getCylinderVector().front()->GetFirstBead()->getID()
        << " " << filament->getCylinderVector().back()->GetSecondBead()->getID() << std::endl;
        
        ///print coordinates
        for (auto cylinder : filament->getCylinderVector()){
            
            auto x = cylinder->GetFirstBead()->force;
            _outputFile<<x[0]<<" "<<x[1]<<" "<<x[2]<<" ";
            
        }
        auto x = filament->getCylinderVector().back()->GetSecondBead()->force;
        _outputFile<<x[0]<<" "<<x[1]<<" "<<x[2];
        _outputFile << std::endl;
    }
    _outputFile <<std::endl;
}

void Output::printStresses(int step) {}

