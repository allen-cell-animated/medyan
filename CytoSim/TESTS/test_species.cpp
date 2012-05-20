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

using namespace std;
using namespace chem;

TEST(SpeciesTest, CTors) {
    //1st CTor
    SpeciesType AType{"Arp2/3",SType::Diffusing};
    Species A{AType, 25};
    EXPECT_EQ(AType,A.getType());
    EXPECT_EQ(25,A.getN());
    
    //2nd CTor
    Species B{"Arp2/3",SType::Diffusing,25};
    EXPECT_EQ(AType,B.getType());
    EXPECT_EQ(25,B.getN());
}

//TEST(SpeciesTest, DeletedFunctions) {
//    Species B{"Arp2/3",SType::Diffusing,25};
//    Species C{"G-Actin",SType::Diffusing,25};
//    C=B; // this should not compile
//    Species D(B); // this should not compile
//}

TEST(SpeciesTest, CopyNumber) {
    Species B{"Arp2/3",SType::Diffusing,25};
    B.setN(32);
    EXPECT_EQ(32,B.getN());
}

TEST(SpeciesTest, SpeciesType) {
    Species B{"Arp2/3",SType::Diffusing,25};
    SpeciesType BType{"Arp2/3",SType::Diffusing};
    EXPECT_EQ(BType,B.getType());
    EXPECT_TRUE(B.is_of_species_type("Arp2/3",SType::Diffusing));
}

TEST(SpeciesTest, Clone) {
    Species B{"Arp2/3",SType::Diffusing,25};
    unique_ptr<Species> C = B.clone();
    EXPECT_EQ(0,C->getN()) << "cloning should initialize copy number to 0";
    EXPECT_EQ(B.getType(),C->getType());
}

TEST(SpeciesTest, Iostream) {
    Species B{"Arp2/3",SType::Diffusing,25};
    ostringstream outStream;
    
    outStream << B;
    string res = outStream.str();
    string b_fullname = B.getFullName() + "[" + std::to_string(B.getN()) + "]";
    EXPECT_EQ(b_fullname,res);
}
