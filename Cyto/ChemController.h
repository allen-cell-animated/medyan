//
//  ChemController.h
//  Cyto
//
//  Created by James Komianos on 7/30/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#ifndef __Cyto__ChemController__
#define __Cyto__ChemController__

#include <iostream>
#include "ChemSimImpl.h"
#include "ChemNRMImpl.h"
#include "ChemGillespieImpl.h"
#include "ChemSimpleGillespieImpl.h"


using namespace chem;

class ChemController {
    
    
    
    void initialize(std::string name) {
        
        ChemSimImpl* csi;
        if(name == "NRM")
            csi = new ChemNRMImpl();
        else if(name == "Gillespie")
            csi = new ChemGillespieImpl();
        else if(name == "SimpleGillespie")
            csi = new ChemSimpleGillespieImpl();
        else {
            std::cout<< "Chem implementation not found. Exiting" <<std::endl;
            exit(EXIT_FAILURE);
        }
        ChemSim::setInstance(csi);
        
        
        ChemSim::initialize();
        
        
        
    }
    
    
    
    
    
    
    
    
    
};


















#endif /* defined(__Cyto__ChemController__) */
