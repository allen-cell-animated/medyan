//
//  example_misc.cpp
//  CytoSim
//
//  Created by Garegin Papoian on 5/29/12.
//  Copyright (c) 2012 University of Maryland. All rights reserved.
//

#include <iostream>

#include <iostream>
#include <fstream>

#include <boost/heap/pairing_heap.hpp>

#include "System.h"
#include "ChemNRMImpl.h"
#include "ChemSim.h"
#include "Signaling.h"


using namespace std;
using namespace chem;



void print_species (RSpecies *s, int i){
    cout << "print_species callback: " << s->getFullName() << ", copy_n=" << s->getN() << endl;
}

struct PrintSpecies {
    void operator() (RSpecies *s, int i){
        cout << "PrintSpecies callback: i=" << _count << "\n";
        cout << (*s);
        ++_count;
    }
    int _count;
};

int main(int argc, const char * argv[])
{
    
    System S;
    
    SpeciesBulk A1("A1",  25);
    SpeciesBulk A2("A2", 25);
    SpeciesBulk A3("A3", 25);
    SpeciesBulk A4("A4", 25);
    SpeciesBulk A5("A5", 25);
    SpeciesBulk A6("A6", 25);
    SpeciesBulk B1("B1", 25);
    
    Reaction r1 = { {&A1,&A3,&A5}, 2, 1, 100.0 };
    Reaction r2 = { {&A2,&A4,&A6}, 2, 1, 50.0 };
    Reaction r3 = { {&A1,&A5,&A6}, 2, 1, 25.0 };
    Reaction r4 = { {&A3,&A4,&A5}, 2, 1, 10.0 };
    Reaction r5 = { {&A1,&A2,&A4}, 2, 1, 5.0 };
    Reaction r6 = { {&B1,&B1,&A1}, 1, 2, 0.05 };
    
    for(int i=0; i<10; ++i)
        r1.makeStep();
    for(int i=0; i<10; ++i)
        r2.makeStep();
    for(int i=0; i<10; ++i)
        r3.makeStep();
    
    cout << "\nTesting r2.getAffectedReactions():\n";
    auto ar = r2.getAffectedReactions();
    cout << "The following reaction," << endl;
    cout << r2 << endl;
    cout << "Affects these reactions:" << endl;
    for(auto r: ar)
        cout << r << endl;
    
    cout << "\n\nTesting computePropensity() " << endl;
    cout << r2 << endl;
    cout << "Current a=" << r2.computePropensity() << endl;
    
    
    
    ChemSignal sm;
    A1.makeSignaling(sm);
    PrintSpecies ps;
    std::function<void (RSpecies *, int)> psf(ps);
    //    sm.connect(A1, print_species);
    //    sm.connect(A1, PrintSpecies());
    //    boost::signals2::shared_connection_block conn_a1(sm.connect(A1,psf), false);
    boost::signals2::shared_connection_block conn_a1(sm.connect(&A1, [](RSpecies *s, int i){cout << *s << endl;}), false);
    
    A2.makeSignaling(sm);
    std::function<void (RSpecies *, int)> psff = [](RSpecies *s, int i){cout << *s << endl;};
    sm.connect(&A2,psff);
    
    
    
    ChemNRMImpl chem_nrm_impl(sm);
    ChemSim chem(&chem_nrm_impl);
    chem.addReaction(&r1);
    chem.addReaction(&r2);
    chem.addReaction(&r3);
    chem.addReaction(&r4);
    chem.addReaction(&r5);
    chem.addReaction(&r6);
    
    // cout << "Num of Reactions is: " << chem.getSize() << ", current time is " << chem.getTime() << endl;
    
    chem.initialize();
    
    cout << "Before starting the NRM runs" << endl;
    chem.printReactions();
    chem.run(100);   
    
    A2.stopSignaling(sm);
    conn_a1.block();
    //    sm.disconnect(A1,psf);
    chem.run(100);   
    
    cout << "Final result:" << endl;
    chem.printReactions();
    
    //    PQNode xx = heap.top();
    //    ReactionNodeNRM *yy = xx.ReactionNodeNRM;
    //
    //    cout << "Pointer sizes, ReactionNodeNRM: " << sizeof (*yy)  << ", PQNode: " << sizeof xx << ", heap size, " << heap.size() << endl;
    //    
    //    yy->printSelf();
    //    yy->makeStep(0.0);
    //    yy->updateHeap();
    //    yy->printSelf();
    //    yy->printDependents();
    
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
    
    cout << "\n" << endl;
    RSpecies &RA1 = A1.getRSpecies();
    cout << "Pointer sizes, Species, RSpecies, Reaction, " 
    << sizeof A1 << ", " <<  sizeof RA1 << ", " << sizeof r1 << "\n\n" << endl;
    
    return 0;
}
