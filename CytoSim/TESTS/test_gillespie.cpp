//
//  test_gillespie.cpp
//  CytoSim
//
//  Created by Garegin Papoian on 7/14/12.
//  Copyright (c) 2012 University of Maryland. All rights reserved.
//

#include <iostream>

#include "gtest/gtest.h"

#include "Species.h"
#include "Reaction.h"
#include "ChemNRMImpl.h"
#include "ChemSim.h"

using namespace std;
using namespace chem;

TEST(ChemNRMTest, SimpleStoichiometryInvariants) {
    SpeciesBulk A1("A1",  100);
    SpeciesBulk A2("A2", 0);
    Reaction r1 = { {&A1,&A2}, 1, 1, 100.0 };

    ChemNRMImpl chem_nrm_impl;
    ChemSim chem(&chem_nrm_impl);
    chem.addReaction(&r1);
    chem.initialize();
//    chem.printReactions();
    chem.run(30);   
//    chem.printReactions();
    EXPECT_EQ(70,A1.getN());
    EXPECT_EQ(30,A2.getN());
}

TEST(ChemNRMTest, StoichiometryInvariants) {
    SpeciesBulk A1("A1",  100);
    SpeciesBulk A2("A2", 0);
    SpeciesBulk A3("A3", 0);
    Reaction r1 = { {&A1,&A2}, 1, 1, 10.0 };
    Reaction r2 = { {&A2,&A1}, 1, 1, 15.0 };
    Reaction r3 = { {&A1,&A3}, 1, 1, 20.0 };
    
    ChemNRMImpl chem_nrm_impl;
    ChemSim chem(&chem_nrm_impl);
    chem.addReaction(&r1);
    chem.addReaction(&r2);
    chem.addReaction(&r3);
    chem.initialize();
    //    chem.printReactions();
    chem.run(30);   
    //    chem.printReactions();
    EXPECT_EQ(100,A1.getN()+A2.getN()+A3.getN());
}

TEST(ChemNRMTest, SimpleSteadyState) {
    SpeciesBulk A1("A1",  100);
    SpeciesBulk A2("A2", 0);
    // A1 <-> A2 with the same forward and backward rates; [A]~[B] at steady state
    Reaction r1 = { {&A1,&A2}, 1, 1, 100.0 };
    Reaction r2 = { {&A2,&A1}, 1, 1, 100.0 };
    
    ChemNRMImpl chem_nrm_impl;
    ChemSim chem(&chem_nrm_impl);
    chem.addReaction(&r1);
    chem.addReaction(&r2);
    chem.initialize();
    //    chem.printReactions();
    chem.run(1000);   
    chem.printReactions();
}