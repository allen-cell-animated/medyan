
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

#include "DRController.h"

#include "Parser.h"
#include "RateChangerImpl.h"

#include "Linker.h"
#include "MotorGhost.h"
#include "Cylinder.h"
#include "BranchingPoint.h"

void DRController::initialize(DynamicRateType& drTypes) {
    
    //filament polymerization changer
    int filamentIndex = 0;
    
    if(SysParams::Chemistry().numFilaments != 0) {
    
        for(auto &changer : drTypes.dFPolymerizationType) {
        
            if(changer == "BROWRATCHET") {
                //get params
                floatingpoint a = SysParams::DynamicRates().dFilPolymerizationCharLength[filamentIndex];
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

    //branching point unbinding changer
    int branchIndex = 0;
    int charLengthIndexbr = 0;
    int charFIndexbr = 0;
    
    if(sum(SysParams::Chemistry().numBrancherSpecies) !=0) {
        
        for(auto &changer : drTypes.dBUnbindingType) {
            
            if(changer == "SLIP") {
                
                //if user did not specify enough parameters, return
                if(charLengthIndexbr >= SysParams::DynamicRates().dBranchUnbindingCharLength.size() )
                    return;
                
                //get the param
                floatingpoint x1 = SysParams::DynamicRates().dBranchUnbindingCharLength[charLengthIndexbr];
                
                //add the rate changer
                BranchingPoint::_unbindingChangers.push_back(new BranchSlip(branchIndex, x1));
                charLengthIndexbr += 1;
            }
            else if(changer == "SLIPF"){
                //if user did not specify enough parameters, return
                if(charFIndexbr >= SysParams::DynamicRates().dBranchUnbindingCharForce
                                                .size() )
                    return;

                //get the param
                floatingpoint x1 = SysParams::DynamicRates()
                        .dBranchUnbindingCharForce[charLengthIndexbr];

                //add the rate changer
                BranchingPoint::_unbindingChangers.push_back(new BranchSlipF(branchIndex, x1));
                charFIndexbr += 1;

            }
            else {
                cout << "Branching point unbinding rate changing form not recognized. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
            
            branchIndex++;
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
                floatingpoint a1 = SysParams::DynamicRates().dLinkerUnbindingAmplitude[ampIndex];
                floatingpoint a2 = SysParams::DynamicRates().dLinkerUnbindingAmplitude[ampIndex + 1];
                
                //now char length
                floatingpoint x1 = SysParams::DynamicRates().dLinkerUnbindingCharLength[charLengthIndex];
                floatingpoint x2 = SysParams::DynamicRates().dLinkerUnbindingCharLength[charLengthIndex + 1];
                
                //add the rate changer
                Linker::_unbindingChangers.push_back(new LinkerCatchSlip(linkerIndex, a1, a2, x1, x2));
                
                charLengthIndex += 2;
                ampIndex += 2;
            }
            
            else if(changer == "SLIP") {
                
                //if user did not specify enough parameters, return
                if(charLengthIndex >= SysParams::DynamicRates().dLinkerUnbindingCharLength.size() )
                    return;
                
                //get the param
                floatingpoint x1 = SysParams::DynamicRates().dLinkerUnbindingCharLength[charLengthIndex];
                
                //add the rate changer
                Linker::_unbindingChangers.push_back(new LinkerSlip(linkerIndex, x1));
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
                floatingpoint f = SysParams::DynamicRates().dMotorUnbindingCharForce[forceIndex];
                
                //add the rate changer
                MotorGhost::_unbindingChangers.push_back(new LowDutyMotorCatch(motorIndex, f));
                forceIndex++;
            }
            else if(changer == "HIGHDUTYCATCH") {
                
                //if user did not specify enough parameters, return
                if(forceIndex >= SysParams::DynamicRates().dMotorUnbindingCharForce.size())
                    return;
                
                //get param
                floatingpoint f = SysParams::DynamicRates().dMotorUnbindingCharForce[forceIndex];
                
                //add the rate changer
                MotorGhost::_unbindingChangers.push_back(new HighDutyMotorCatch(motorIndex, f));
                forceIndex++;
            }
            else if(changer == "MOTORSLIP") {
                
                //if user did not specify enough parameters, return
                if(forceIndex >= SysParams::DynamicRates().dMotorUnbindingCharForce.size())
                    return;
                
                //get param
                floatingpoint f = SysParams::DynamicRates().dMotorUnbindingCharForce[forceIndex];
                
                //add the rate changer
                MotorGhost::_unbindingChangers.push_back(new MotorSlip(motorIndex, f));
                forceIndex++;
            }

            else if(changer == "LOWDUTYCATCHSLIP") {
                cout << "Catch-slip bond implementation of low duty motor not complete. Exiting." << endl;
                exit(EXIT_FAILURE);
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
                floatingpoint f = SysParams::DynamicRates().dMotorWalkingCharForce[forceIndex];
                
                //add the rate changer
                MotorGhost::_walkingChangers.push_back(new LowDutyMotorStall(motorIndex, 0, f));
                forceIndex++;
            }
            else if(changer == "HIGHDUTYSTALL") {
                
                //if user did not specify enough parameters, return
                if(forceIndex >= SysParams::DynamicRates().dMotorWalkingCharForce.size())
                    return;
                
                //get the param
                floatingpoint f = SysParams::DynamicRates().dMotorWalkingCharForce[forceIndex];
                
                //add the rate changer
                MotorGhost::_walkingChangers.push_back(new HighDutyMotorStall(motorIndex, 0, f));
                forceIndex++;
            }
            else if(changer == "TWOHEADSTALL") {
                 
                 //if user did not specify enough parameters, return
                 if(forceIndex >= SysParams::DynamicRates().dMotorWalkingCharForce.size())
                     return;
                 
                 //get the param
                 floatingpoint f = SysParams::DynamicRates().dMotorWalkingCharForce[forceIndex];
                 
                 float walkingrate = SysParams::DynamicRates().dMotorWalkingRate[forceIndex];
                 
                 //add the rate changer
                 MotorGhost::_walkingChangers.push_back(new TwoHeadStall(motorIndex, 0, f, walkingrate));
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

