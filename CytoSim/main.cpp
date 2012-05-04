//
//  main.cpp
//  CytoSim-Experimenting
//
//  Created by Garegin Papoian on 4/19/12.
//  Copyright (c) 2012 University of Maryland. All rights reserved.
//

#include <iostream>
#include "System.h"
#include "ChemNRM.h"
#include <boost/heap/pairing_heap.hpp>

using namespace std;


int main(int argc, const char * argv[])
{
    Species Arp23("Arp2/3", SType::Diffusing, 25);
    Species Actin("Actin", SType::Diffusing, 30);
    Arp23.printSelf();
    Actin.printSelf();
//    std::array<Species*,1> lsp = {&Arp23};
//    std::array<Species*,1> rsp = {&Actin};
    //Reaction<1,1> r1{1.1,2.2,lsp,rsp};
//    Reaction<1,1> r1{1.14,2.222,Arp23,Actin};
//    r1.setFRate(3.33);
//    r1.setBRate(4.44);
//    r1.printSelf();
//    r1.makeFStep();
//    r1.printSelf();
//    cout << "Pointer sizes, Species vs Reaction<1,1>" << sizeof Arp23 << " " << sizeof r1 << endl;
//    
//    ReactionNode rn = {r1,true};
//    rn.setTau(23.4);
//    vector<ReactionNode> vrn = {rn};
//    vector<ReactionNode> vrn2;
//    //vrn2.push_back(rn);
//    
//    boost::heap::pairing_heap<ReactionNode> pq;
//    
//    cout << "Pointer sizes, ReactionNode: " << sizeof vrn[0] << endl;
//    cout << rn.getReactBase() << " " << vrn[0].getReactBase() << " " << vrn2[0].getReactBase() << endl;
//    rn.makeStep();
//    vrn[0].makeStep();
//    vrn2[0].makeStep();
//    r1.printSelf();
//    cout << "Propensities: " << r1.getFPropensity() << ", " << r1.getBPropensity() << endl;
    return 0;
}

