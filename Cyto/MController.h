//
//  MController.h
//  Cyto
//
//  Created by James Komianos on 8/4/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#ifndef __Cyto__MController__
#define __Cyto__MController__

#include "SubSystem.h"
#include "Mcommon.h"
#include <iostream>


class MController {
    
public:
    
    void run(System* ps, std::string solver) {
     
        if(solver == "FletcherRieves")
            FletcherRievesMethod(ps);
        else if (solver == "PolakRibiere")
            PolakRibiereMethod(ps);
        
        else {
            std::cout<< "Mechanical algorithm not found. Exiting." <<std::endl;
            exit(EXIT_FAILURE);
        }
    }
    
    
    
};




#endif /* defined(__Cyto__MController__) */
