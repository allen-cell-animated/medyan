
//------------------------------------------------------------------
//  **M3SYM** - Simulation Package for the Mechanochemical
//              Dynamics of Active Networks, 3rd Generation
//
//  Copyright (2014) Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the Papoian lab page for installation and documentation:
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

#include "MathFunctions.h"

using namespace mathfunc;

void BasicSnapshot::print(int step) {
    
    _outputFile.precision(10);
    
    //print first line (step number, time, number of filaments, linkers, motors, branchers)
    _outputFile << step << " " << tau() << " " <<
        FilamentDB::instance()->size() << " " <<
        LinkerDB::instance()->size() << " " <<
        MotorGhostDB::instance()->size() << " " <<
        BranchingPointDB::instance()->size() << endl;
    
    for(auto &filament : *FilamentDB::instance()) {
        
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
    
    
    for(auto &linker : *LinkerDB::instance()) {
        
        //print first line
        _outputFile << "L " << linker->getLinkerID()<< " " <<
            linker->getLinkerType() << endl;
        
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

    for(auto &motor : *MotorGhostDB::instance()) {
        
        //print first line
        _outputFile << "M " << motor->getMotorID() << " " <<
            motor->getMotorType() << endl;
        
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
    
    for(auto &branch : *BranchingPointDB::instance()) {
        
        //print first line
        _outputFile << "B " << branch->getBranchID() << " " <<
        branch->getBranchType() << endl;
        
        //print coordinates
        auto x = branch->coordinate;
        _outputFile<<x[0]<<" "<<x[1]<<" "<<x[2] << endl;
    }
    
    _outputFile <<endl;
}

void BirthTimes::print(int step) {
    
    _outputFile.precision(10);
    
    //print first line (step number, time, number of filaments, linkers, motors, branchers)
    _outputFile << step << " " << tau() << " " <<
        FilamentDB::instance()->size() << " " <<
        LinkerDB::instance()->size() << " " <<
        MotorGhostDB::instance()->size() << " " <<
        BranchingPointDB::instance()->size() << endl;
    
    for(auto &filament : *FilamentDB::instance()) {
        
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
        _outputFile<< filament->getCylinderVector().back()->getSecondBead()->getBirthTime();
        _outputFile << endl;
    }
    for(auto &linker : *LinkerDB::instance()) {
        
        //print first line
        _outputFile << "L " << linker->getLinkerID()<< " " <<
        linker->getLinkerType() << endl;
        
        //print birth times
        _outputFile << linker->getBirthTime() << endl;
    }
    
    for(auto &motor : *MotorGhostDB::instance()) {
        
        //print first line
        _outputFile << "M " << motor->getMotorID() << " " <<
        motor->getMotorType() << endl;
        
        //print birth times
        _outputFile << motor->getBirthTime() << endl;
    }
    
    for(auto &branch : *BranchingPointDB::instance()) {
        
        //print first line
        _outputFile << "B " << branch->getBranchID() << " " <<
        branch->getBranchType() << endl;
        
        //print birth times
        _outputFile << branch->getBirthTime() << endl;
    }
    
    _outputFile <<endl;
}

void Forces::print(int step) {
    
    _outputFile.precision(10);
    
    //print first line (step number, time, number of filaments, linkers, motors, branchers)
    _outputFile << step << " " << tau() << " " <<
        FilamentDB::instance()->size() << " " <<
        LinkerDB::instance()->size() << " " <<
        MotorGhostDB::instance()->size() << " " <<
        BranchingPointDB::instance()->size() << endl;
    
    for(auto &filament : *FilamentDB::instance()) {
        
        //print first line(Filament ID, length, left_delta, right_delta
        _outputFile << "F " << filament->getID() << " " <<
        filament->getCylinderVector().size() + 1 << " " <<
        filament->getDeltaMinusEnd() << " " << filament->getDeltaPlusEnd() << endl;
        
        //print force
        for (auto cylinder : filament->getCylinderVector()){
            
            double forceMag= cylinder->getFirstBead()->calcForceSquare();
            forceMag = sqrt(forceMag);
            _outputFile<<forceMag << " ";
            
        }
        //print last bead force
        double forceMag = filament->getCylinderVector().back()->
                          getSecondBead()->calcForceSquare();
        forceMag = sqrt(forceMag);
        _outputFile<<forceMag;
        
        _outputFile << endl;
    }
    
    for(auto &linker : *LinkerDB::instance()) {
        
        //print first line
        _outputFile << "L " << linker->getLinkerID()<< " " <<
        linker->getLinkerType() << endl;
        
        //print stretch force
        _outputFile << 0.0 << endl;
    }
    
    for(auto &motor : *MotorGhostDB::instance()) {
        
        //print first line
        _outputFile << "M " << motor->getMotorID() << " " <<
        motor->getMotorType() << endl;
        
        //print stretch force
        _outputFile << 0.0 << endl;
    }
    
    for(auto &branch : *BranchingPointDB::instance()) {
        
        //print first line
        _outputFile << "B " << branch->getBranchID() << " " <<
        branch->getBranchType() << endl;
        
        //print bending force
        _outputFile << 0.0 << endl;
    }
    
    _outputFile <<endl;
}


void Stresses::print(int step) {}

