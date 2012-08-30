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

namespace chem {

    class Visitor;

class Composite : public Component {
private:
    std::vector<std::unique_ptr<Component>> _children;

public:
    Composite() :  Component() {}
    virtual ~Composite() noexcept {}
    
    virtual bool apply (Visitor &v);
    virtual bool apply_if (ConditionalVisitor &v);
    
    virtual bool isComposite() {return true;}
    
    virtual bool numberOfChildren() {return children().size();}

    virtual std::string getFullName() const {return "Composite";}; 
    
    virtual void addChild (std::unique_ptr<Component> &&child) {
        _children.push_back(std::move(child));
        _children.back()->setParent(this);
    }
    
    virtual void removeChild (Component* child) {
        auto child_iter = std::find_if(_children.begin(),_children.end(),
                    [&child](const std::unique_ptr<Component> &element)
                     {
                         return element.get()==child ? true : false;
                     });
    }
    
    virtual std::vector<std::unique_ptr<Component>>& children () {return _children;}

    virtual const std::vector<std::unique_ptr<Component>>& children () const {return _children;}

    virtual Component* children (size_t i) {return _children[i].get();}
    
};
       
} // end of chem
#endif
