
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

void DRController::initialize(DynamicRateType& drTypes) {
    
    //filament polymerization changer
    int filamentIndex = 0;
    
    if(SysParams::Chemistry().numFilaments != 0) {
    
        for(auto &changer : drTypes.dFPolymerizationType) {
        
            if(changer == "BROWRATCHET") {
                //get params
                double a = SysParams::DynamicRates().dFilPolymerizationCharLength[filamentIndex];
                Cylinder::_polyChanger.push_back(new BrownianRatchet(a));
            }
            else if(changer == "") {}
            else {
                cout << "Filament polymerization rate changing form not recognized. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
            filamentIndex++;
        }
    }

    //linker unbinding changer
    int charLengthIndex = 0;
    int ampIndex = 0;
    int linkerIndex = 0;
    
    if(sum(SysParams::Chemistry().numLinkerSpecies) != 0) {
    
        for(auto &changer : drTypes.dLUnbindingType) {
            
            if(changer == "CATCHSLIP") {
                
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
                Linker::_unbindingChangers.push_back(new CatchSlip(linkerIndex, a1, a2, x1, x2));
                
                charLengthIndex += 2;
                ampIndex += 2;
            }
            
            else if(changer == "SLIP") {
                
                //if user did not specify enough parameters, return
                if(charLengthIndex >= SysParams::DynamicRates().dLinkerUnbindingCharLength.size() )
                    return;
                
                //get the param
                double x1 = SysParams::DynamicRates().dLinkerUnbindingCharLength[charLengthIndex];
                
                //add the rate changer
                Linker::_unbindingChangers.push_back(new Slip(linkerIndex, x1));
                charLengthIndex += 1;
            }
            else {
                cout << "Linker unbinding rate changing form not recognized. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
            
            linkerIndex++;
        }
        
    }
    int forceIndex = 0;
    int motorIndex = 0;
    
    if(sum(SysParams::Chemistry().numMotorSpecies) != 0) {
    
        //motor unbinding changer
        for(auto &changer : drTypes.dMUnbindingType) {
            
            if(changer == "LOWDUTYCATCH") {
                
                //if user did not specify enough parameters, return
                if(forceIndex >= SysParams::DynamicRates().dMotorUnbindingCharForce.size())
                    return;
                
                //get param
                double f = SysParams::DynamicRates().dMotorUnbindingCharForce[forceIndex];
                
                //add the rate changer
                MotorGhost::_unbindingChangers.push_back(new LowDutyCatch(motorIndex, f));
                forceIndex++;
            }
            if(changer == "LOWDUTYCATCHSLIP") {
                
                //if user did not specify enough parameters, return
                if(forceIndex >= SysParams::DynamicRates().dMotorUnbindingCharForce.size())
                    return;
                
                //get param
                double fCatch = SysParams::DynamicRates().dMotorUnbindingCharForce[forceIndex];
                double fSlip  = SysParams::DynamicRates().dMotorUnbindingCharForce[forceIndex + 1];
                
                //add the rate changer
                MotorGhost::_unbindingChangers.push_back(new LowDutyCatchSlip(motorIndex, fCatch, fSlip));
                forceIndex += 2;
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

            if(changer == "LOWDUTYSTALL") {
                
                //if user did not specify enough parameters, return
                if(forceIndex >= SysParams::DynamicRates().dMotorWalkingCharForce.size())
                    return;
                
                //get the param
                double f = SysParams::DynamicRates().dMotorWalkingCharForce[forceIndex];
                
                //add the rate changer
                MotorGhost::_walkingChangers.push_back(new LowDutyStall(motorIndex, 0, f));
                forceIndex++;
            }
            else {
                cout << "Motor walking rate changing form not recognized. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
            motorIndex++;
        }
    }
}

