
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

#include "DRController.h"

#include "Parser.h"
#include "RateChangerImpl.h"

#include "Linker.h"
#include "MotorGhost.h"
#include "Cylinder.h"

void DRController::initialize(DynamicRateTypes& drTypes) {
    
    //filament polymerization changer
    if(drTypes.dFPolymerizationType == "BROWRATCHET") {
        //get params
        double a = SysParams::DynamicRates().dFilPolymerizationCharLength;
        Cylinder::_polyChanger = new BrownianRatchet(a);
    }
    else if(drTypes.dFPolymerizationType == "") {}
    else {
        cout << "Filament polymerization rate changing form not recognized. Exiting." << endl;
        exit(EXIT_FAILURE);
    }
    
    //linker unbinding changer
    int charLengthIndex = 0;
    int ampIndex = 0;
    int linkerIndex = 0;
    
    for(auto &changer : drTypes.dLUnbindingType) {
        
        if(changer == "BASICCATCHSLIP") {
            
            //if user did not specify enough parameters, return
            if(ampIndex + 1 >= SysParams::DynamicRates().dLinkerUnbindingAmplitude.size() ||
               charLengthIndex + 1 >= SysParams::DynamicRates().dLinkerUnbindingCharLength.size())
                return;
            
            //get two params for each, amp
            double a1 = SysParams::DynamicRates().dLinkerUnbindingAmplitude[ampIndex];
            double a2 = SysParams::DynamicRates().dLinkerUnbindingAmplitude[ampIndex + 1];
            
            //now char length
            double x1 = SysParams::DynamicRates().dLinkerUnbindingCharLength[charLengthIndex];
            double x2 = SysParams::DynamicRates().dLinkerUnbindingCharLength[charLengthIndex + 1];
            
            //add the rate changer
            Linker::_unbindingChangers.push_back(
                new BasicCatchSlip(linkerIndex, a1, a2, x1, x2));
            
            charLengthIndex += 2;
            ampIndex += 2;
        }
        
        else if(changer == "BASICSLIP") {
            
            //if user did not specify enough parameters, return
            if(charLengthIndex >= SysParams::DynamicRates().dLinkerUnbindingCharLength.size() )
                return;
            
            //get the param
            double x1 = SysParams::DynamicRates().dLinkerUnbindingCharLength[charLengthIndex];
            
            //add the rate changer
            Linker::_unbindingChangers.push_back(new BasicSlip(linkerIndex, x1));
            charLengthIndex += 1;
        }
        else {
            cout << "Linker unbinding rate changing form not recognized. Exiting." << endl;
            exit(EXIT_FAILURE);
        }
        
        linkerIndex++;
    }
    int forceIndex = 0;
    int motorIndex = 0;
    
    //motor unbinding changer
    for(auto &changer : drTypes.dMUnbindingType) {
        
        if(changer == "LOWDUTYPCMCATCH") {
            
            //if user did not specify enough parameters, return
            if(forceIndex >= SysParams::DynamicRates().dMotorUnbindingCharForce.size())
                return;
            
            //get param
            double f = SysParams::DynamicRates().dMotorUnbindingCharForce[forceIndex];
            
            //add the rate changer
            MotorGhost::_unbindingChangers.push_back(new LowDutyPCMCatch(motorIndex, f));
            forceIndex++;
        }
        else {
            cout << "Motor unbinding rate changing form not recognized. Exiting." << endl;
            exit(EXIT_FAILURE);
        }
        motorIndex++;
    }
    forceIndex = 0;
    motorIndex = 0;
    
    //motor walking 
    for(auto &changer : drTypes.dMWalkingType) {

        if(changer == "LOWDUTYHILLSTALL") {
            
            //if user did not specify enough parameters, return
            if(forceIndex >= SysParams::DynamicRates().dMotorWalkingCharForce.size())
                return;
            
            //get the param
            double f = SysParams::DynamicRates().dMotorWalkingCharForce[forceIndex];
            
            //add the rate changer
            MotorGhost::_walkingChangers.push_back(new LowDutyHillStall(motorIndex, f));
            forceIndex++;
        }
        else {
            cout << "Motor walking rate changing form not recognized. Exiting." << endl;
            exit(EXIT_FAILURE);
        }
        motorIndex++;
    }
}

