
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
#include "RateChanger.h"

#include "Linker.h"
#include "MotorGhost.h"
#include "Cylinder.h"

void DRController::initialize(DynamicRateTypes& drTypes) {
    
    //filament polymerization changer
    if(drTypes.dFPolymerizationType == "BROWRATCHET") {
        //get params
        double a = SystemParameters::DynamicRates().dFilPolymerizationCharLength;
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
    
    for(auto &changer : drTypes.dLUnbindingType) {
        
        if(changer == "CATCHSLIP") {
            
            //if user did not specify enough parameters, return
            if(ampIndex + 1 >= SystemParameters::DynamicRates().
                               dLinkerUnbindingAmplitude.size() ||
               charLengthIndex + 1 >= SystemParameters::DynamicRates().
                                      dLinkerUnbindingCharLength.size() )
                return;
            
            //get two params for each
            double a1 = SystemParameters::DynamicRates().
                        dLinkerUnbindingAmplitude[ampIndex];
            double a2 = SystemParameters::DynamicRates().
                        dLinkerUnbindingAmplitude[ampIndex + 1];
            double x1 = SystemParameters::DynamicRates().
                        dLinkerUnbindingCharLength[charLengthIndex];
            double x2 = SystemParameters::DynamicRates().
                        dLinkerUnbindingCharLength[charLengthIndex + 1];
            
            //add the rate changer
            Linker::_unbindingChangers.push_back(new CatchSlipBond(a1, a2, x1, x2));
            charLengthIndex += 2;
            ampIndex += 2;
        }
        
        else if(changer == "SLIP") {
            
            //if user did not specify enough parameters, return
            if(charLengthIndex >= SystemParameters::DynamicRates().
                                  dLinkerUnbindingCharLength.size() )
                return;
            
            //get the param
            double x1 = SystemParameters::DynamicRates().
                        dLinkerUnbindingCharLength[charLengthIndex];
            
            //add the rate changer
            Linker::_unbindingChangers.push_back(new SlipBond(x1));
            charLengthIndex += 1;
        }
        else {
            cout << "Linker unbinding rate changing form not recognized. Exiting." << endl;
            exit(EXIT_FAILURE);
        }
    }
    charLengthIndex = 0;
    ampIndex = 0;
    
    //motor unbinding changer
    for(auto &changer : drTypes.dMUnbindingType) {
        
        if(changer == "CATCHSLIP") {
            
            //if user did not specify enough parameters, return
            if(ampIndex + 1 >= SystemParameters::DynamicRates().
                               dMotorUnbindingAmplitude.size() ||
               charLengthIndex + 1 >= SystemParameters::DynamicRates().
                               dMotorUnbindingCharLength.size() )
                return;
            
            //get two params for each
            double a1 = SystemParameters::DynamicRates().
                        dMotorUnbindingAmplitude[ampIndex];
            double a2 = SystemParameters::DynamicRates().
                        dMotorUnbindingAmplitude[ampIndex + 1];
            double x1 = SystemParameters::DynamicRates().
                        dMotorUnbindingCharLength[charLengthIndex];
            double x2 = SystemParameters::DynamicRates().
                        dMotorUnbindingCharLength[charLengthIndex + 1];
            
            //add the rate changer
            MotorGhost::_unbindingChangers.push_back(new CatchSlipBond(a1, a2, x1, x2));
            charLengthIndex += 2;
            ampIndex += 2;
        }
        
        else if(changer == "SLIP") {
            
            //if user did not specify enough parameters, return
            if(charLengthIndex >= SystemParameters::DynamicRates().
                                  dMotorUnbindingCharLength.size() )
                return;
            
            //get the param
            double x1 = SystemParameters::DynamicRates().
                        dMotorUnbindingCharLength[charLengthIndex];
            
            //add the rate changer
            MotorGhost::_unbindingChangers.push_back(new SlipBond(x1));
            charLengthIndex += 1;
        }
        else {
            
            cout << "Motor unbinding rate changing form not recognized. Exiting." << endl;
            exit(EXIT_FAILURE);
        }
    }
    
    charLengthIndex = 0;
    
    //motor walking 
    for(auto &changer : drTypes.dMWalkingType) {

        if(changer == "EXPSTALL") {
            
            //if user did not specify enough parameters, return
            if(charLengthIndex >= SystemParameters::DynamicRates().
                                  dMotorWalkingCharLength.size() )
                return;
            
            //get the param
            double x1 = SystemParameters::DynamicRates().
                        dMotorWalkingCharLength[charLengthIndex];
            
            //add the rate changer
            MotorGhost::_walkingChangers.push_back(new ExpStall(x1));
            charLengthIndex += 1;
        }
        else {
            cout << "Motor walking rate changing form not recognized. Exiting." << endl;
            exit(EXIT_FAILURE);
        }
    }
}

