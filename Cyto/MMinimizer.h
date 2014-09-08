//
//  MMinimizer.h
//  Cyto
//
//  Created by Konstantin Popov on 9/8/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#ifndef __Cyto__MMinimizer__
#define __Cyto__MMinimizer__

#include <iostream>

class MController;

class Minimizer
{

public:
    
    void virtual Equlibrate(MController*) = 0;
    
};




#endif /* defined(__Cyto__MMinimizer__) */
