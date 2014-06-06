//
//  test_compartment.cpp
//  CytoSim
//
//  Created by James Komianos on 6/5/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//


#define DO_THIS_TEST

#ifdef DO_THIS_TEST

#include <iostream>
#include "gtest/gtest.h"
#include "Compartment.h"
#include "CompartmentContainer.h"
#include "ReactionBase.h"

using namespace std;
using namespace chem;


//COMPARTMENT TESTS
//Testing neighbor capabilities

TEST(CompartmentTest, Neighbors)
{
    Compartment *C1 = new Compartment();
    Compartment *C2 = C1->clone();
    
    C1->addNeighbour(C2);
    C2->addNeighbour(C1);
    
    EXPECT_EQ(1U, C1->numberOfNeighbours());
    EXPECT_EQ(1U, C2->numberOfNeighbours());
    
    Compartment *C3 = C2->clone();
    
    C1->addNeighbour(C3);
    C3->addNeighbour(C1);
    
    EXPECT_EQ(2U, C1->numberOfNeighbours());
    EXPECT_EQ(1U, C3->numberOfNeighbours());
    
    C1->removeNeighbour(C2);
    EXPECT_EQ(1U, C1->numberOfNeighbours());
}

//Testing diffusion rates

TEST(CompartmentTest, DiffusionRates)
{
    Compartment *C1 = new Compartment();
    Species *actin = C1->addSpecies("Actin",99U);
    C1->setDiffusionRate(actin,2000);
    
    Compartment *C2 = C1->clone();
    
    C1->addNeighbour(C2);
    C2->addNeighbour(C1);
    
    C1->generateDiffusionReactions();
    C2->generateDiffusionReactions();
    
    EXPECT_EQ(C1->getDiffusionRate("Actin"), 2000);
    EXPECT_EQ(C2->getDiffusionRate("Actin"), 2000);
    
    C1->setDiffusionRate(actin, 9000);
    EXPECT_EQ(C1->getDiffusionRate("Actin"), 9000);
    
}

//COMPARTMENTGRID TESTS

//Testing initialize
TEST(CompartmentContainerTest, initialize)
{
    CompartmentGrid<3> g{50,50,50};
    CompartmentSpatial<3> &Cproto = g.getProtoCompartment();
    Species *M1 = Cproto.addSpecies("Myosin",1U);
    Cproto.setDiffusionRate(M1,2000);
    Species *M2 = Cproto.addSpecies("Fascin",6U);
    Cproto.setDiffusionRate(M2,2000);
    Cproto.addInternal<Reaction,1,1>({M1,M2}, 40.2);
    Cproto.addInternal<Reaction,1,1>({M2,M1}, 90.9);
    
    vector<float> sides{100.0,100.0,100.0};
    Cproto.setSides(sides.begin());
    g.initialize();
    
    EXPECT_EQ(g.countSpecies(), 250000);
    EXPECT_EQ(g.countReactions(), 1720000);
    
}

//Testing neighbors
TEST(CompartmentContainerTest, Neighbors)
{
    CompartmentGrid<3> g{50,50,50};
    CompartmentSpatial<3> &Cproto = g.getProtoCompartment();
    Species *M1 = Cproto.addSpecies("Myosin",1U);
    Cproto.setDiffusionRate(M1,2000);
    Species *M2 = Cproto.addSpecies("Fascin",6U);
    Cproto.setDiffusionRate(M2,2000);
    
    vector<float> sides{100.0,100.0,100.0};
    Cproto.setSides(sides.begin());
    
    Cproto.addInternal<Reaction,1,1>({M1,M2}, 40.2);
    Cproto.addInternal<Reaction,1,1>({M2,M1}, 90.9);
    g.initialize();
    
    CompartmentSpatial<3>* C1 = g.getCompartment(10U,10U,10U);
    CompartmentSpatial<3>* C2 = g.getCompartment(10U,11U,10U);
    
    EXPECT_EQ(C1->coords()[0], 10*100.0);
    EXPECT_EQ(C1->coords()[1], 10*100.0);
    EXPECT_EQ(C1->coords()[2], 10*100.0);
    
    EXPECT_EQ(C1->numberOfNeighbours(), 6U);
    
    C1->removeNeighbour(C2);
    EXPECT_EQ(C1->numberOfNeighbours(), 5U);
    C1->addNeighbour(C2);
}


#endif //DO_THIS_TEST