//
//  VisitorConcrete.h
//  CytoSim
//
//  Created by Garegin Papoian on 6/2/12.
//  Copyright (c) 2012 University of Maryland. All rights reserved.
//

#ifndef CytoSim_VisitorConcrete_h
#define CytoSim_VisitorConcrete_h

#include "Visitor.h"
#include "Composite.h"
#include "utility.h"

namespace chem {

/// ConcreteVisitor is an example visitor, following the regular Visitor pattern
class ConcreteVisitor : public Visitor {
protected:
    virtual bool visitImpl(Component *c) override {
        std::cout << "Visiting, " << c->getFullName() << ", which has " 
        << c->countSpecies() << ", Species overall. " << std::endl;
        return true;
    }
};

///// FindFirstSpeciesVisitor is an example visitor, which is propagated through the
///// Composite hierarchy until the first Composite node is found which contains Species
///// with a certain name
//class FindFirstSpeciesVisitor : public Visitor {
//private: 
//    std::string _name;
//public:
//    FindFirstSpeciesVisitor(const std::string &name) : _name(name) {}
//protected:
//    virtual bool visitImpl(Component *c) override {
//        Composite *C = dynamic_cast<Composite*>(c);
//        for(auto &s : C->species()){
//            std::cout << "Visiting Species:" << "\n" << (*s) << std::endl;
//            if(s->getName() == _name){
//                std::cout << "Found Species " << (*s) << std::endl;
//                return false;
//            }
//        }
//        return true;
//    }
//};
    
class CompositeVisitor : public Visitor {
protected:
    virtual bool visitImpl(Component *c) override {
        std::cout << "Visiting, " << c->getFullName() << ", which has " 
        << c->countSpecies() << ", Species overall. " << std::endl;
        return true;
    }
    virtual bool predImpl(Component *c) override {
        return isSame<Composite>(*c);
    }

};
    
//class ProtoFilamentVisitor : public ConditionalVisitor {
//public:
//    virtual bool visit(Component *c) {
//        std::cout << "Visiting, " << c->getFullName() << ", which has " 
//        << c->countSpecies() << ", Species overall. " << std::boolalpha << isSame<ProtoFilament>(*c) << std::endl;
//        return true;
//    }
//    virtual bool pred(Component *c) {
//        return typeid(*c) == typeid(ProtoFilament);
//    }
//    
//};

} // end of chem

#endif
