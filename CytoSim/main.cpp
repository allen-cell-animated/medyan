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

Bead make_bead(Species *s){
    Bead bead(s);
    bead.x()=3.0;
    bead.y()=5.0;
    bead.z()=7.0;    
    return bead;
}

int main(int argc, const char * argv[])
{
    System S;

    S.parseSpeciesProtos();
    S.parseSystem();

    Species *arp23_proto = S.SpeciesProto("Arp2/3",SType::Diffusing);
    Species *gactin_proto = S.SpeciesProto("G-Actin",SType::Diffusing);
    Species *motor_proto = S.SpeciesProto("Motor",SType::Diffusing);
    Species *profilin_proto = S.SpeciesProto("Profilin",SType::Diffusing);
        
    Species r1(arp23_proto->getType()); r1.setN(20);
    Species p1(gactin_proto->getType()); p1.setN(30);

    std::array<Species*,2> sp2{{&r1,&p1}};
//    Reaction<1,1> r1p1 {12.0, 5.0, sp2};
//    
//    for(int i=0;i<6;++i)
//        r1p1.doStep(true);
//    r1p1.printSelf();

    Species r2(profilin_proto->getType()); r2.setN(40);
    Species p2(motor_proto->getType()); p2.setN(50);
    std::array<Species*,4> sp4 = {{&r1,&r2,&p1,&p2}};
//    Reaction<2,2> r2p2{0,0,sp4};
//    r2p2.printSelf();
//    for(int i=0;i<7;++i)
//        r2p2.doStep(true);
//    r2p2.printSelf();

//    cout << "main(x..): sizeof(r1): " << sizeof(r1) << endl;
//    cout << "main(x.x): sizeof(r1p1): " << sizeof(r1p1) << endl;

    Bead bb1(&r1);
    bb1.x()=2.2;
    bb1.printSelf();

    bb1 = make_bead(&p1);    
    bb1.printSelf();
    
    Bead bb2{&r2,{{1.4,2.4,3.4}}};
    bb2.printSelf();
    
    cout << "Still alive in main..." << endl;
    
//    S.parseSystem();
//    Space1D* s1d = new Space1D();
//    S.setSpace(s1d);
//    S.setSpaceOptions(std::make_tuple(100,0,0));
//    S.initializeCompartments();    
//    cout << S.getNumCompartments() << " " << sizeof(int) << endl;

    return 0;
}

