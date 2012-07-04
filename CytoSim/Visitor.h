//
//  Visitor.h
//  CytoSim
//
//  Created by Garegin Papoian on 6/2/12.
//  Copyright (c) 2012 University of Maryland. All rights reserved.
//

#ifndef CytoSim_Visitor_h
#define CytoSim_Visitor_h

#include <iostream>
#include "Component.h"

namespace chem {

class Visitor {
public:
    virtual bool visit(Component *c) = 0;
    virtual ~Visitor() {}
};
    
class ConditionalVisitor {
public:
    virtual bool visit_if(Component *c) {
        std::cout << "ConditionalVisitor::visit_if: checking " << c->getFullName() << std::endl;
        if(pred(c)){
            return visit(c);
        }
        else{
            std::cout << "ConditionalVisitor::visit_if: Predicate failed, moving on..." << std::endl;
            return true;
        }
    }
    virtual bool visit(Component *c) = 0;
    virtual bool pred(Component *c) = 0;
    virtual ~ConditionalVisitor() {}
};
    
    
} // end of chem
#endif
