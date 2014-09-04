//
//  MComponent.h
//  CytoMech
//
//  Created by Konstantin Popov on 6/16/14.
//  Copyright (c) 2014 Konstantin Popov. All rights reserved.
//

#ifndef __CytoMech__MComponent__
#define __CytoMech__MComponent__

#include <iostream>

class MComposite;

class MComponent{
    
private:
    MComposite *_parent; ///< The parent of this object. May be a nullptr if this object has no parent.
    
public:
    /// Default Constructor; Parent is assigned to nullptr
    MComponent() : _parent(nullptr) {}
    
    /// Virtual Destructor
    virtual ~MComponent() {};
    /// Returns the pointer to the parent node. The returned value could be a nullptr if a parent does not exist.
    virtual MComposite* getParent() {return _parent;};
};

#endif /* defined(__CytoMech__MComponent__) */
