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
    Species A1("A1", SType::Diffusing, 25);
    Species A2("A2", SType::Diffusing, 26);
    Species A3("A3", SType::Diffusing, 27);
    Species A4("A4", SType::Diffusing, 28);
    Species A5("A5", SType::Diffusing, 29);

    Species B1("B1", SType::Diffusing, 30);
    Species B2("B2", SType::Diffusing, 31);
    Species B3("B3", SType::Diffusing, 32);
    Species B4("B4", SType::Diffusing, 33);
    Species B5("B5", SType::Diffusing, 34);
    
    Species C1("C1", SType::Diffusing, 35);
    Species C2("C2", SType::Diffusing, 36);
    Species C3("C3", SType::Diffusing, 37);
    Species C4("C4", SType::Diffusing, 38);
    Species C5("C5", SType::Diffusing, 39);
    
    Reaction<2,1> r1 = { {&A1,&B1,&C1}, 3.0 };
    Reaction<2,1> r2 = { {&A1,&B2,&C2}, 3.0 };
    Reaction<2,1> r3 = { {&A1,&B3,&C3}, 3.0 };
    Reaction<2,1> r4 = { {&A1,&B4,&C4}, 3.0 };
    Reaction<2,1> r5 = { {&A1,&B5,&C5}, 3.0 };

    Reaction<2,1> r6 = { {&B2,&C2,&A1}, 3.0 };

    
    r1.printSelf();
    cout << endl;
    
    for(int i=0; i<10; ++i)
        r1.makeStep();
//    r1.makeStep();
    r1.printSelf();
    cout << "Currtent a=" << r1.computePropensity() << endl;
    cout << "Pointer sizes, Species vs Reaction<1,1>" << sizeof A1 << " " << sizeof r1 << "\n\n" << endl;
    
    for(auto r = A1.beginBReactions(); r!=A1.endBReactions(); ++r)
        (*r)->printSelf();
    
    cout << "\nTesting r1.getAffectedReactions():\n";
    auto ar = r6.getAffectedReactions();
    for(auto r: ar)
        r->printSelf();
    
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

