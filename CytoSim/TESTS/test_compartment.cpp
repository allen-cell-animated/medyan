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

using namespace std;
using namespace chem;








// Testing the creation of diffusion reactions between compartments
// Testing the compartment neighbors
// Testing the get compartment functionality
TEST(CompartmentContainerTest, Main) {
    
    CompartmentsSimpleGrid<3> g{50,50,50};
    Compartment &CProto = g.getProtoCompartment();
    Species *M1 = Cproto.addSpecies("Myosin",1U);
    Cproto.setDiffusionRate(M1,2000);
    Species *M2 = Cproto.addSpecies("Fascin",6U);
    Cproto.setDiffusionRate(M2,2000);
    ReactionBase *RM1M2 = Cproto.addInternal<Reaction,1,1>({M1,M2}, 40.2);
    ReactionBase *RM2M1 = Cproto.addInternal<Reaction,1,1>({M2,M1}, 90.9);
    g.initialize();
    
    EXPECT_EQ(g.countSpecies(), 128);
    EXPECT_EQ(g.countReactions(), 704);
    
    
    
    
    
    
    
}






























#endif //DO_THIS_TEST