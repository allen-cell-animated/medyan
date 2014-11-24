//
//  Component.cpp
//  CytoSim
//
//  Created by Garegin Papoian on 5/31/12.
//  Copyright (c) 2012 University of Maryland. All rights reserved.
//
#include <iostream>

#include "Component.h"

#include "Composite.h"
#include "Visitor.h"

Composite* Component::getRoot() {
    if(isRoot())
        return static_cast<Composite*>(this);
    else
        return getParent()->getRoot();
}

bool Component::apply (Visitor &v) {
    return v.visit(this);
}