//
//  Component.h
//  CytoSim
//
//  Created by Garegin Papoian on 5/31/12.
//  Copyright (c) 2012 University of Maryland. All rights reserved.
//

#ifndef CytoSim_Component_h
#define CytoSim_Component_h

namespace chem {
    
    class Composite;
    class Visitor;
    class ConditionalVisitor;
    
    class Component {
    private:
        Composite *_parent;
    public:
        Component() : _parent(nullptr) {}
        Composite* getParent() {return _parent;}
        void setParent (Composite *other) {_parent=other;}
        bool isRoot() const {return _parent==nullptr? true : false;}
        Composite* getRoot();
        virtual bool hasChildren() {return false;}
        virtual ~Component() noexcept {}
        virtual std::string getFullName() const = 0;
        virtual size_t countSpecies() const = 0;
        virtual bool apply (Visitor &v);
        virtual bool apply_if (ConditionalVisitor &v);
    };

}// end of chem


#endif
