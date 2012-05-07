//
//  main.cpp
//  CytoSim-Experimenting
//
//  Created by Garegin Papoian on 4/19/12.
//  Copyright (c) 2012 University of Maryland. All rights reserved.
//

#include <iostream>
#include <boost/heap/pairing_heap.hpp>

#include "System.h"
#include "ChemNRMImpl.h"

using namespace std;

int main(int argc, const char * argv[])
{
    SpeciesContainer sc;
    sc.addSpecies("A1", "A1", SType::Diffusing, 25);
    sc.addSpecies("A2", "A2", SType::Diffusing, 25);
    sc.addSpecies("A3", "A3", SType::Diffusing, 25);
    sc.addSpecies("A4", "A4", SType::Bulk, 25);
    sc.addSpecies("A5", "A5", SType::Bulk, 25);
    sc.addSpecies("A6", "A6", SType::Bulk, 25);
    
    Species* A1 = sc.getSpecies("A1");
    Species* A2 = sc.getSpecies("A2");
    Species* A3 = sc.getSpecies("A3");
    Species* A4 = sc.getSpecies("A4");
    Species* A5 = sc.getSpecies("A5");
    Species* A6 = sc.getSpecies("A6");
            
    Reaction r1 = { {A1,A3,A5}, 2, 1, 3.2 };
    Reaction r2 = { {A2,A4,A6}, 2, 1, 3.5 };
    Reaction r3 = { {A1,A5,A6}, 2, 1, 3.7 };
    Reaction r4 = { {A3,A4,A5}, 2, 1, 3.9 };
    Reaction r5 = { {A1,A2,A4}, 2, 1, 4.1 };
    
    for(int i=0; i<10; ++i)
        r1.makeStep();
    for(int i=0; i<10; ++i)
        r2.makeStep();
    for(int i=0; i<10; ++i)
        r3.makeStep();
    
    cout << "\nTesting r2.getAffectedReactions():\n";
    auto ar = r2.getAffectedReactions();
    cout << "The following reaction," << endl;
    r2.printSelf();
    cout << "Affects these reactions:" << endl;
    for(auto r: ar)
        r->printSelf();
    
    cout << "\n\nTesting computePropensity() " << endl;
    r2.printSelf();
    cout << "Current a=" << r2.computePropensity() << endl;
    cout << "Pointer sizes, Species vs Reaction, " << sizeof (*A1) << " " << sizeof r1 << "\n\n" << endl;
    
    for(auto r = A2->beginReactantReactions(); r!=A2->endReactantReactions(); ++r)
        (*r)->printSelf();
    

    boost_heap heap;
    RNode rn1(r1,heap);
    RNode rn2(r2,heap);
    RNode rn3(r3,heap);
    // Testing destructors
    {
        RNode rn4(r4,heap);
        RNode rn5(r5,heap);
        rn5.printSelf();
        rn5.printDependents();
    }
    

    PQNode xx = heap.top();
    RNode *yy = xx.rnode;

    cout << "Pointer sizes, ReactionNodeNRM: " << sizeof (*yy)  << ", PQNode: " << sizeof xx << ", heap size, " << heap.size() << endl;
    
    yy->printSelf();
    yy->makeStep(0.0);
    yy->updateHeap();
    yy->printSelf();
    yy->printDependents();
    
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

