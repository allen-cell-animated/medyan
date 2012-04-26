//
//  main.cpp
//  CytoSim-Experimenting
//
//  Created by Garegin Papoian on 4/19/12.
//  Copyright (c) 2012 University of Maryland. All rights reserved.
//

#include <iostream>
#include "System.h"


using std::cout;
using std::endl;

int main(int argc, const char * argv[])
{
    System S;

    S.parseSpeciesTypes();
    S.parseSystem();

    SpeciesType *arp23 = S.getSpeciesType("Arp2/3");
    SpeciesType *gactin = S.getSpeciesType("G-Actin");
    SpeciesType *profilin = S.getSpeciesType("Profilin");
    SpeciesType *motor = S.getSpeciesType("Motor");
    
    Species r1(*arp23); r1.setN(20);
    Species p1(*gactin); p1.setN(30);
    std::array<Species*,2> sp2{{&r1,&p1}};
    Reaction<1,1> r1p1 {12.0, 5.0, sp2};
    for(int i=0;i<6;++i)
        r1p1.doStep(true);
    r1p1.printSelf();
    
    Species r2(*profilin); r2.setN(40);
    Species p2(*motor); p2.setN(50);
    std::array<Species*,4> sp4 = {{&r1,&r2,&p1,&p2}};
    Reaction<2,2> r2p2{0,0,sp4};
    r2p2.printSelf();
    for(int i=0;i<7;++i)
        r2p2.doStep(true);
    r2p2.printSelf();
    
    cout << "main(x..): sizeof(r1): " << sizeof(r1) << endl;
    cout << "main(x.x): sizeof(r1p1): " << sizeof(r1p1) << endl;
//    S.parseSystem();
//    Space1D* s1d = new Space1D();
//    S.setSpace(s1d);
//    S.setSpaceOptions(std::make_tuple(100,0,0));
//    S.initializeCompartments();    
//    cout << S.getNumCompartments() << " " << sizeof(int) << endl;
    
    return 0;
}

