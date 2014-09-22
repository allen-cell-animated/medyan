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

#include "ChemInitializerImpl.h"
#include "CFilamentElement.h"
#include "common.h"

///SimpleInitializer is a concrete implementation of the ChemInitailizerImpl class, which sets up CCylinders
///to have simple actin network interactions
class SimpleInitializerImpl : public ChemInitializerImpl {

private:
    ///REACTION RATES
    //basic
    float _k_on_plus = 21.0;
    float _k_on_minus = 0.0;
    float _k_off_plus = 1.4;
    float _k_off_minus = 1.4;
    
    //Diffusion rate
    float _diffusion_rate = 2000.0;
    
    
public:
    ///Initialize the compartment grid
    virtual void initializeGrid();

    ///Initializer
    ///@param length - starting length of the CCylinder initialized
    ///@param species - list of species to initialize in CCylinder
    virtual CCylinder* createCCylinder(Filament* pf, Compartment* c, bool extensionFront, bool extensionBack);
    
    ///Remove a CCylinder
    virtual void removeCCylinder(Filament* pf, bool retractionFront, bool retractionBack);
    
};

#endif /* defined(__Cyto__SimpleInitializerImpl__) */
