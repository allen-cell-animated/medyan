//
//  SimpleInitializerImpl.h
//  Cyto
//
//  Created by James Komianos on 7/30/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#ifndef __Cyto__SimpleInitializerImpl__
#define __Cyto__SimpleInitializerImpl__

#include <iostream>

#include "CompartmentContainer.h"
#include "CCylinder.h"
#include "CFilamentElement.h"
#include "Mcommon.h"
#include "common.h"


namespace chem {

    ///SimpleInitializer is a concrete implementation of the ChemInitailizerImpl class, which sets up CCylinders
    ///to have simple actin network interactions
    class SimpleInitializerImpl : public ChemInitializerImpl {

    private:
        ///REACTION RATES
        //basic
        float _k_on_plus = 21.0;
        float _k_off_plus = 1.4;
        float _k_off_minus = 1.4;
        
        //capping
        float _k_capping_on_plus = 50.0;
        float _k_capping_off_plus = 0.06;
        
        //formin
        float _k_formin_on_plus = 10.0;
        float _k_formin_off_plus = 1.4;
        float _k_accel_on_plus = 100.0;
        
        //Diffusion rate
        float _diffusion_rate = 2000.0;
        
        
    public:
        ///Initializer, based on the given simulation
        ///@param length - starting length of the CCylinder initialized
        ///@param species - list of species to initialize in CCylinder
        virtual CCylinder* createCCylinder(Compartment* c, std::vector<std::string> species, int length);
        
        ///Remove a CCylinder, based on the given simulation
        virtual void removeCCylinder(CCylinder *cylinder);
        
    };

}; // end namespace chem


#endif /* defined(__Cyto__SimpleInitializerImpl__) */
