//
//  MComposite.h
//  CytoMech
//
//  Created by Konstantin Popov on 6/16/14.
//  Copyright (c) 2014 Konstantin Popov. All rights reserved.
//

#ifndef __CytoMech__MComposite__
#define __CytoMech__MComposite__

#include <iostream>
#include "MComponent.h"

class Bead;

class MComposite : public MComponent{
public:
    /// Default Constructor does nothing.
    MComposite() :  MComponent() {}
    
    /// Virtual Destructor does nothing.
    virtual ~MComposite() noexcept {}
    virtual void DeleteBead(Bead* pb) = 0;
};

#endif /* defined(__CytoMech__MComposite__) */
