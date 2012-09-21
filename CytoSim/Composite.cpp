//
//  Composite.cpp
//  CytoSim
//
//  Created by Garegin Papoian on 5/29/12.
//  Copyright (c) 2012 University of Maryland. All rights reserved.
//

#include <iostream>
#include "Composite.h"
#include "Visitor.h"

using namespace chem;
    

bool Composite::apply (Visitor &v) {
    bool res_self = v.visit(this); //pre-order
    if(!res_self)
        return false;
    for (auto &c : children()) {
        bool res_child = c->apply(v);
        if(!res_child)
            return false;
    }
    return true;
}


bool Composite::apply (SpeciesVisitor &v) {
    bool res_self = apply_impl(v); //pre-order
    if(!res_self)
        return false;
    for (auto &c : children()) {
        bool res_child = c->apply(v);
        if(!res_child)
            return false;
    }
    return true;
}

bool Composite::apply (ReactionVisitor &v) {
    bool res_self = apply_impl(v); //pre-order
    if(!res_self)
        return false;
    for (auto &c : children()) {
        bool res_child = c->apply(v);
        if(!res_child)
            return false;
    }
    return true;
}