//
//  Reactable.h
//  Cyto
//
//  Created by James Komianos on 11/25/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#ifndef __Cyto__Reactable__
#define __Cyto__Reactable__

#include <iostream>

///Reactable class is for a reactable object in the subsystem
class Reactable {
    
public:
    virtual void updateReactionRates() = 0;
};

#endif /* defined(__Cyto__Reactable__) */
