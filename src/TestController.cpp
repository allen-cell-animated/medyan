#ifdef CAMKII

//
//  TestController.cpp
//  MEDYAN
//
//  Created by James Komianos on 7/27/17.
//  Copyright © 2017 University of Maryland. All rights reserved.
//

#include "TestController.h"
#include "SubSystem.h"

#include "common.h"
#include "Camkii.h"

#include "Output.h"
#include "MController.h"
#include "GController.h"
#include "CController.h"
#include "DRController.h"

void TestController::run() {
    
#ifdef CHEMISTRY
    double tauLastSnapshot = 0;
    double tauLastMinimization = 0;
    double tauLastNeighborList = 0;
    double oldTau = 0;
    
    long stepsLastSnapshot = 0;
    long stepsLastMinimization = 0;
    long stepsLastNeighborList = 0;
    
    long totalSteps = 0;
#endif

#ifdef CHEMISTRY
     _subSystem->updateBindingManagers();
#endif
#ifdef DYNAMICRATES
        updateReactionRates();
#endif
    cout << "TEST CONTROLLER" << endl;
    cout << "Current simulation time = "<< tau() << endl;
    //restart phase ends
#ifdef CHEMISTRY
    tauLastSnapshot = tau();
    oldTau = 0;
#endif
    for(auto o: _outputs) o->print(0);
    
    cout << "Starting simulation..." << endl;
    
    int i = 1;

    //MAKE NEW CAMKII
    //PERFORM MECHANICAL MINIMIZATION
    vector<double> coord = {0,0,0};
    //Bead* next = _subSystem->addTrackable<Bead>(coord, this, 1);
    Camkii* c = _subSystem->addTrackable<Camkii>(_subSystem, 0, coord);
    
    c->printSelf();
    
    _mController->run();
    updatePositions();
    
    
    //if runtime was specified, use this
    if(!areEqual(_runTime, 0.0)) {
        
#ifdef CHEMISTRY
        while(tau() <= _runTime) {
            //run ccontroller
            if(!_cController->run(_minimizationTime)) {
                for(auto o: _outputs) o->print(i);
                break;
            }
            
            //add the last step
            tauLastSnapshot += tau() - oldTau;
            tauLastMinimization += tau() - oldTau;
            tauLastNeighborList += tau() - oldTau;
#endif
#if defined(MECHANICS) && defined(CHEMISTRY)
            //run mcontroller, update system
            if(tauLastMinimization >= _minimizationTime) {
                _mController->run();
                updatePositions();
                
                tauLastMinimization = 0.0;
            }
            
            if(tauLastSnapshot >= _snapshotTime) {
                cout << "Current simulation time = "<< tau() << endl;
                for(auto o: _outputs) o->print(i);
                i++;
                tauLastSnapshot = 0.0;
            }
#elif defined(MECHANICS)
            for(auto o: _outputs) o->print(i);
            i++;
#endif
            
#ifdef DYNAMICRATES
            updateReactionRates();
#endif
            
#ifdef CHEMISTRY
            // update neighbor lists
            if(tauLastNeighborList >= _neighborListTime) {
                updateNeighborLists();
                tauLastNeighborList = 0.0;
            }
            
            //move the boundary
            moveBoundary(tau() - oldTau);
            
            //special protocols
            executeSpecialProtocols();
            
            oldTau = tau();
        }
#endif
    }
    //if run steps were specified, use this
    if(_runSteps != 0) {
        
#ifdef CHEMISTRY
        while(totalSteps <= _runSteps) {
            //run ccontroller
            if(!_cController->runSteps(_minimizationSteps)) {
                for(auto o: _outputs) o->print(i);
                break;
            }
            
            //add the last step
            stepsLastSnapshot += _minimizationSteps;
            stepsLastMinimization += _minimizationSteps;
            stepsLastNeighborList += _minimizationSteps;
            
            totalSteps += _minimizationSteps;
#endif
#if defined(MECHANICS) && defined(CHEMISTRY)
            //run mcontroller, update system
            if(stepsLastMinimization >= _minimizationSteps) {
                _mController->run();
                updatePositions();
                
                stepsLastMinimization = 0;
            }
            
            if(stepsLastSnapshot >= _snapshotSteps) {
                cout << "Current simulation time = "<< tau() << endl;
                for(auto o: _outputs) o->print(i);
                i++;
                stepsLastSnapshot = 0;
            }
#elif defined(MECHANICS)
            for(auto o: _outputs) o->print(i);
            i++;
#endif
            
#ifdef DYNAMICRATES
            updateReactionRates();
#endif
            
#ifdef CHEMISTRY
            // update neighbor lists
            if(stepsLastNeighborList >= _neighborListSteps) {
                updateNeighborLists();
                stepsLastNeighborList = 0;
            }
            
            //move the boundary
            moveBoundary(tau() - oldTau);
            
            //special protocols
            executeSpecialProtocols();
        }
#endif
    }
    
    //print last snapshots
    for(auto o: _outputs) o->print(i);
    
    cout << "Done with simulation!" << endl;
}
#endif