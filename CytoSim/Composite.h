//
//  Composite.h
//  CytoSim
//
//  Created by Garegin Papoian on 5/29/12.
//  Copyright (c) 2012 University of Maryland. All rights reserved.
//

#ifndef CytoSim_Composite_h
#define CytoSim_Composite_h

#include <string>
#include <vector>
#include <deque>
#include <typeinfo>

#include "Component.h"
#include "Species.h"
#include "Visitor.h"

namespace chem {
    
class Composite : public Component {
private:
    std::vector<std::unique_ptr<Composite>> _children;

public:
    Composite() :  Component() {}
    
    virtual bool apply (Visitor &v) {
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
    
    virtual bool apply_if (ConditionalVisitor &v) {
        bool res_self = v.visit_if(this); //pre-order
        if(!res_self) 
            return false;
        for (auto &c : children()) {
            bool res_child = c->apply_if(v);
            if(!res_child) 
                return false;
        }
        return true;
    }

    virtual ~Composite() noexcept {}
    
    virtual bool hasChildren() {return children().size()>0 ? true : false;}

    virtual std::string getFullName() const {return "Composite";}; 
    
    virtual void addChild (std::unique_ptr<Composite> &&child) {
        _children.push_back(std::move(child));
        _children.back()->setParent(this);
    }
    
    virtual std::vector<std::unique_ptr<Composite>>& children () {return _children;}

    virtual const std::vector<std::unique_ptr<Composite>>& children () const {return _children;}

    virtual Composite* children (size_t i) {return _children[i].get();}
    
};
       
} // end of chem
#endif
