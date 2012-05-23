//
//  test_species.cpp
//  CytoSim
//
//  Created by Garegin Papoian on 5/17/12.
//  Copyright (c) 2012 University of Maryland. All rights reserved.
//

// Note: This test omits many functions of Species that interact with Reaction objects. 
//        Separate tests will cover those methods.

#include <iostream>
#include "gtest/gtest.h"
#include "Species.h"
#include "System.h"

using namespace std;
using namespace chem;

TEST(SpeciesTest, CTors) {
    System S;
    Species *A = S.addSpecies("Arp2/3",SType::Bulk,25);
    
    //1st CTor
    SpeciesType AType{"Arp2/3",SType::Bulk};
    EXPECT_EQ(AType,A->getType());
    EXPECT_EQ(25,A->getN());
}

//TEST(SpeciesTest, DeletedFunctions) {
//    Species B{"Arp2/3",SType::Bulk,25};
//    Species C{"G-Actin",SType::Bulk,25};
//    C=B; // this should not compile
//    Species D(B); // this should not compile
//}

TEST(SpeciesTest, CopyNumber) {
    System S;
    Species *B = S.addSpecies("Arp2/3",SType::Bulk,25);
    B->setN(32);
    EXPECT_EQ(32,B->getN());
}

TEST(SpeciesTest, SpeciesType) {
    System S;
    Species *B = S.addSpecies("Arp2/3",SType::Bulk,25);

    SpeciesType BType{"Arp2/3",SType::Bulk};
    EXPECT_EQ(BType,B->getType());
    EXPECT_TRUE(B->is_of_species_type("Arp2/3",SType::Bulk));
}

TEST(SpeciesTest, Clone) {
    System S;
    SpeciesBulk *B = S.addSpecies("Arp2/3",SType::Bulk,25);

    System S2;
    Species* C = B->clone(S2);
    EXPECT_EQ(B->getN(),C->getN()) << "cloning should initialize copy number to 0";
    EXPECT_EQ(B->getType(),C->getType());
}

TEST(SpeciesTest, Iostream) {    
    System S;
    Species *B = S.addSpecies("Arp2/3",SType::Bulk,25);

    ostringstream outStream;
    
    outStream << *B;
    string res = outStream.str();
    string b_fullname = B->getFullName() + "[" + std::to_string(B->getN()) + "]";
    EXPECT_EQ(b_fullname,res);
}
