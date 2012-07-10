//
//  example_composite.cpp
//  CytoSim
//
//  Created by Garegin Papoian on 7/10/12.
//  Copyright (c) 2012 University of Maryland. All rights reserved.
//

#include <iostream>
#include <fstream>

#include <boost/heap/pairing_heap.hpp>

#include "System.h"
#include "ChemNRMImpl.h"
#include "ChemSim.h"
#include "Signaling.h"

#include "Composite.h"
#include "VisitorConcrete.h"

using namespace std;
using namespace chem;

int main(int argc, const char * argv[])
{
    auto C = make_unique<Composite>();
    cout << C->getFullName() << endl;
    C->pushBackSpeciesUnique(make_unique<SpeciesBulk>("G-Actin",42));
    C->pushBackSpecies<SpeciesBulk>("Profilin",31);
    
    auto a1 = C->species(0);
    auto a2 = C->species(1);
    
    cout << *a1 << "\n" << *a2 << endl;
    
    for(auto &s: C->species()){
        cout << s->getFullName() << s->getN() << endl; 
        cout << &s << endl;
    }
    
    cout << "\n\n\n";
    
    auto D = make_unique<Composite>();
    D->pushBackSpecies<SpeciesBulk>("ATP",3300);
    C->addChild(std::move(D));
    
    Composite *F = C->getRoot();
    F->pushBackSpecies<SpeciesBulk>("Arp2/3",22);
    
    
//    auto pf = make_unique<ProtoFilament>("FA"); //F-Actin
//    pf->growPlusEnd(10);
//    C->addChild(std::move(pf));
    
    for(auto &s: F->species()){
        cout << s->getFullName() << s->getN() << endl; 
        cout << &s << endl;
    }
    
    cout << "F has " << F->countSpecies() << " species" << endl;
    
    cout << "F is of class Component? " << boolalpha << isSame<Component>(*F) << endl; 
    cout << "F is of class Composite? " << boolalpha << isSame<Composite>(*F) << endl; 
    cout << "F is of class Composite? " << boolalpha << isSame<Species>(*F) << endl; 
    
    ConcreteVisitor cv;
    F->apply(cv);
    
    FindFirstSpeciesVisitor ffsv("Arp2/3");
    F->apply(ffsv);
    
    CompositeVisitor comp_vis;
    F->apply_if(comp_vis);
    //    ConcreteConditionalVisitor ccv([](Component *c){return c->hasChildren() ? true : false;});
    
    cout << "\n";
    
//    ProtoFilamentVisitor pf_vis;
//    F->apply_if(pf_vis);
    
    
    cout << "Main exited..." << endl;
    return 0;
}

