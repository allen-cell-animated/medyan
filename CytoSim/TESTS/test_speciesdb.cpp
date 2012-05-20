//
//  test_speciesdb.cpp
//  CytoSim
//
//  Created by Garegin Papoian on 5/19/12.
//  Copyright (c) 2012 University of Maryland. All rights reserved.
//

#include <iostream>
#include "gtest/gtest.h"
#include "SpeciesDB.h"

using namespace std;
using namespace chem;

class SpeciesDBTest : public ::testing::Test {
protected:
    virtual void SetUp() {
        sdb1.create("A1", "A1", SType::Diffusing, 25);
        sdb1.create("A2", "A2", SType::Bulk, 32);
    }
    // virtual void TearDown() {} // not necessary - everything should self-destruct
    
    SpeciesDB sdb0; // this should be initially empty
    SpeciesDB sdb1; // this should initially have two elements
};

TEST_F(SpeciesDBTest, Size) {
    EXPECT_EQ(0, sdb0.size());
    EXPECT_EQ(2, sdb1.size());
}

TEST_F(SpeciesDBTest, Get) {
    Species *a1 = sdb1.get("A1");
    Species *a2 = sdb1.get("A2");
    EXPECT_EQ("A1{Diffusing}", a1->getFullName());
    EXPECT_EQ(25, a1->getN());
    EXPECT_EQ("A2{Bulk}", a2->getFullName());
    EXPECT_EQ(32, a2->getN());
    
    // Now let's give a bogus key - an exception should be thrown
    EXPECT_ANY_THROW(sdb1.get("A3"));
}

TEST_F(SpeciesDBTest, Create) {
    Species *c1 = sdb1.create("C1", "C1", SType::Walking, 25);
    EXPECT_EQ("C1{Walking}", c1->getFullName());
    EXPECT_EQ(25, c1->getN());

    // Now let's give an existing key - an exception must be thrown
    EXPECT_ANY_THROW(sdb1.create("C1", "X", SType::Diffusing, 42));
}

TEST_F(SpeciesDBTest, Clone) {
    Species *a1 = sdb1.get("A1");
    Species *c1 = sdb1.clone("C1", *a1);
    EXPECT_EQ("A1{Diffusing}", c1->getFullName());
    EXPECT_EQ(0, c1->getN()) << "cloned Species should be initialized to 0 copy number";
}