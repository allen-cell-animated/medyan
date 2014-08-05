//
//  test_compartments.cpp
//  CytoSim
//
//  Created by James Komianos on 6/9/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

//#define DO_THIS_TEST

#ifdef DO_THIS_TEST

#include <iostream>
#include "gtest/gtest.h"
#include "CompartmentContainer.h"

using namespace std;

//Testing neighbor capabilities
TEST(CompartmentTest, Neighbors) {

    Compartment *C1 = new Compartment;
    Compartment *C2 = C1->clone();
    C1->addNeighbour(C2);
    C2->addNeighbour(C1);
    
    EXPECT_EQ(1, C1->numberOfNeighbours());
    EXPECT_EQ(1, C2->numberOfNeighbours());
    
    Compartment *C3 = C2->clone();
    C1->addNeighbour(C3);
    
    EXPECT_EQ(2, C1->numberOfNeighbours());
    EXPECT_EQ(0, C3->numberOfNeighbours());
    
    C1->removeNeighbour(C2);
    
    EXPECT_EQ(1, C1->numberOfNeighbours());
}

//Testing species and reaction generation
TEST(CompartmentTest, SpeciesAndReactions){
    
    Compartment *C1 = new Compartment;
    Species *actin = C1->addSpecies("Actin",99U);
    C1->setDiffusionRate(actin,2000);
    Species *profilin = C1->addSpecies("Profilin",29U);
    C1->setDiffusionRate(profilin,2000);
    Species *arp23 = C1->addSpecies("Arp2/3",33U);
    C1->setDiffusionRate(arp23,2000);
    
    Compartment *C2 = C1->clone();
    
    EXPECT_EQ(2000,C1->getDiffusionRate("Actin"));
    EXPECT_EQ(2000,C2->getDiffusionRate("Actin"));
    
    C1->setDiffusionRate(actin, 9000);
    EXPECT_EQ(9000, C1->getDiffusionRate("Actin"));
    EXPECT_EQ(2000, C2->getDiffusionRate("Actin"));
    
    
    C1->addNeighbour(C2);
    C2->addNeighbour(C1);
    C1->generateAllDiffusionReactions();
    C2->generateAllDiffusionReactions();
    
    EXPECT_EQ(3, C1->numberOfSpecies());
    EXPECT_EQ(3, C2->numberOfSpecies());
    
    EXPECT_EQ(0, C1->numberOfInternalReactions());
    EXPECT_EQ(3, C1->numberOfReactions());
    EXPECT_EQ(0, C2->numberOfInternalReactions());
    EXPECT_EQ(3, C2->numberOfReactions());

}

//Testing neighbors, species and reaciton generation
TEST(CompartmentContainerTest, Main) {
    
    CompartmentGrid<3> ccv{50, 50, 50};
    CompartmentSpatial<3> &Cproto = ccv.getProtoCompartment();
    Species *M1 = Cproto.addSpecies("Myosin",1U);
    Cproto.setDiffusionRate(M1,2000);
    Species *M2 = Cproto.addSpecies("Fascin",6U);
    Cproto.setDiffusionRate(M2,2000);
    vector<float> sides{100.0,100.0,100.0};
    Cproto.setSides(sides.begin());
    
    Cproto.addInternal<Reaction,1,1>({M2,M1}, 90.9);
    Cproto.addInternal<Reaction,1,1>({M1,M2}, 40.2);

    ccv.initialize();
    
    EXPECT_EQ(ccv.countSpecies(), 250000);
    EXPECT_EQ(ccv.countReactions(), 1720000);
    
    CompartmentSpatial<3> *C1 = ccv.getCompartment(10U,10U,10U);
    
    EXPECT_EQ(6, C1->numberOfNeighbours());
    EXPECT_EQ(2, C1->numberOfSpecies());
    EXPECT_EQ(14, C1->numberOfReactions());
    
    CompartmentSpatial<3> *C2 = ccv.getCompartment(10U,11U,10U);
    C1->removeNeighbour(C2);
    
    EXPECT_EQ(5, C1->numberOfNeighbours());
    
    EXPECT_EQ(1000.0, C1->coords()[0]);
    EXPECT_EQ(1100.0, C2->coords()[1]);
    C1->addNeighbour(C2);

}

#endif // DO_THIS_TEST






