#ifdef CAMKII
//
//  TestController.hpp
//  MEDYAN
//
//  Created by James Komianos on 7/27/17.
//  Copyright © 2017 University of Maryland. All rights reserved.
//

#ifndef TestController_hpp
#define TestController_hpp

#include <stdio.h>
#include "Controller.h"


class TestController: public Controller {
    
public:
    
    TestController(SubSystem* s): Controller(s) {};
    ~TestController() {};
    
    
    virtual void run();
};

#endif /* TestController_hpp */
#endif